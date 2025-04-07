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
    Detects if the synthesis involves conversion of a nitrile to a hydroxylamine derivative.
    """
    # Track if we found nitrile to hydroxylamine conversion
    found_nitrile_conversion = False

    def dfs_traverse(node):
        nonlocal found_nitrile_conversion

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if reactants contain nitrile
            reactants = reactants_part.split(".")
            nitrile_present = False
            nitrile_reactant = None

            for reactant in reactants:
                if checker.check_fg("Nitrile", reactant):
                    nitrile_present = True
                    nitrile_reactant = reactant
                    print(f"Found nitrile in reactant: {reactant}")

            # Check if product contains hydroxylamine derivative
            if nitrile_present:
                # Check for specific reaction types first
                if checker.check_reaction(
                    "N-hydroxyimidamide from nitrile and hydroxylamine", rsmi
                ) or checker.check_reaction(
                    "Amidoxime from nitrile and hydroxylamine", rsmi
                ):
                    print(f"Found nitrile to hydroxylamine conversion reaction: {rsmi}")
                    found_nitrile_conversion = True
                else:
                    # Fallback to checking product structure
                    product_mol = Chem.MolFromSmiles(product_part)
                    if product_mol:
                        # Check for N-hydroxyimidamide structure
                        if product_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[C](=N)[N][O]")
                        ) or checker.check_fg("Oxime", product_part):
                            print(
                                f"Found hydroxylamine derivative in product: {product_part}"
                            )
                            # Verify this is a conversion from nitrile to hydroxylamine derivative
                            # by checking that other reactants contain hydroxylamine or related compounds
                            for reactant in reactants:
                                if reactant != nitrile_reactant:
                                    reactant_mol = Chem.MolFromSmiles(reactant)
                                    if reactant_mol and reactant_mol.HasSubstructMatch(
                                        Chem.MolFromSmarts("[N][O]")
                                    ):
                                        print(
                                            f"Found hydroxylamine-like structure in reactant: {reactant}"
                                        )
                                        found_nitrile_conversion = True
                                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return found_nitrile_conversion
