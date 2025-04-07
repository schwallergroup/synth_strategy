#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    Detects if the synthesis involves diazonium chemistry for aromatic functionalization.

    Diazonium chemistry typically involves:
    1. Formation of diazonium salt from primary aromatic amine and nitrite
    2. Subsequent reactions of the diazonium salt (substitution, coupling, etc.)
    """
    has_diazonium_chemistry = False

    # Track molecules that might be involved in diazonium chemistry
    potential_diazonium_precursors = set()

    def dfs_traverse(node, depth=0):
        nonlocal has_diazonium_chemistry

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if this molecule contains a diazonium group
            if "N#N" in mol_smiles or "N=[N+]" in mol_smiles or "[N+]#N" in mol_smiles:
                print(f"Found potential diazonium intermediate: {mol_smiles}")
                has_diazonium_chemistry = True
                return

            # Track aromatic amines as potential diazonium precursors
            if checker.check_fg("Aniline", mol_smiles):
                print(
                    f"Tracking aromatic amine as potential diazonium precursor: {mol_smiles}"
                )
                potential_diazonium_precursors.add(mol_smiles)

        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a known diazo addition reaction
            if checker.check_reaction("Diazo addition", rsmi):
                print(f"Detected diazo addition reaction: {rsmi}")
                has_diazonium_chemistry = True
                return

            # Check for evidence of diazonium chemistry
            has_nitrite_source = False
            has_aromatic_amine = False

            for reactant in reactants:
                # Check for nitrite sources (more specific than just nitro compounds)
                if (
                    "ONO" in reactant
                    or "NO2" in reactant
                    or checker.check_fg("Nitroso", reactant)
                ):
                    has_nitrite_source = True
                    print(f"Found nitrite source: {reactant}")

                # Check for primary aromatic amines
                if checker.check_fg("Aniline", reactant):
                    has_aromatic_amine = True
                    print(f"Found aromatic amine: {reactant}")
                    potential_diazonium_precursors.add(reactant)

            # Check for diazonium formation conditions
            if has_aromatic_amine and has_nitrite_source:
                print(
                    f"Detected potential diazonium formation conditions in reaction: {rsmi}"
                )
                has_diazonium_chemistry = True
                return

            # Check for diazonium intermediates in product
            if "N#N" in product or "N=[N+]" in product or "[N+]#N" in product:
                print(f"Found diazonium group in product: {product}")
                has_diazonium_chemistry = True
                return

            # Check for products of diazonium substitution reactions (Sandmeyer-type)
            if has_aromatic_amine or any(
                precursor in reactants for precursor in potential_diazonium_precursors
            ):
                # Common diazonium substitution products
                if (
                    checker.check_fg("Aromatic halide", product)
                    or checker.check_fg("Phenol", product)
                    or checker.check_fg("Nitrile", product)
                    or checker.check_fg("Azide", product)
                ):
                    print(f"Found potential diazonium substitution product: {product}")
                    has_diazonium_chemistry = True
                    return

                # Check for azo coupling products
                if "-N=N-" in product or "c-N=N-c" in product:
                    print(f"Found potential azo coupling product: {product}")
                    has_diazonium_chemistry = True
                    return

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_diazonium_chemistry
