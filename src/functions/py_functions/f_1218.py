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
    """Check for amide formation from acid chloride in the synthesis route"""
    found = False

    def dfs(node, depth=0):
        nonlocal found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for amide formation reactions
            if (
                checker.check_reaction(
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    rxn_smiles,
                )
                or checker.check_reaction(
                    "Acyl chloride with secondary amine to amide", rxn_smiles
                )
                or checker.check_reaction("Schotten-Baumann_amide", rxn_smiles)
                or checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    rxn_smiles,
                )
                or checker.check_reaction("Acylation of primary amines", rxn_smiles)
                or checker.check_reaction("Acylation of secondary amines", rxn_smiles)
            ):

                # Verify that amide is formed
                try:
                    reactants = rxn_smiles.split(">")[0].split(".")
                    product = rxn_smiles.split(">")[-1]

                    has_acyl_halide = any(
                        checker.check_fg("Acyl halide", r) for r in reactants
                    )
                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Aniline", r)
                        for r in reactants
                    )
                    has_amide = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )

                    if (
                        (has_acyl_halide or any("C(=O)Cl" in r for r in reactants))
                        and (has_amine or any("NH2" in r for r in reactants))
                        and (has_amide or "NC(=O)" in product)
                    ):
                        found = True
                        print(f"Found amide formation at depth {depth}: {rxn_smiles}")
                except Exception as e:
                    print(f"Error checking amide formation: {e}")

            # Additional check for amide formation from carboxylic acid
            if not found:
                try:
                    reactants = rxn_smiles.split(">")[0].split(".")
                    product = rxn_smiles.split(">")[-1]

                    has_carboxylic = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    )
                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Aniline", r)
                        for r in reactants
                    )
                    has_amide = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )

                    if has_carboxylic and has_amine and has_amide:
                        found = True
                        print(
                            f"Found amide formation from carboxylic acid at depth {depth}: {rxn_smiles}"
                        )
                except Exception as e:
                    print(f"Error checking amide formation from carboxylic acid: {e}")

        # Recursively check children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    return found
