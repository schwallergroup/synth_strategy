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
    Detects if the synthesis route involves a urea disconnection strategy.
    """
    urea_disconnection_found = False

    def dfs_traverse(node):
        nonlocal urea_disconnection_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a urea synthesis reaction directly
            urea_reactions = [
                "Urea synthesis via isocyanate and primary amine",
                "Urea synthesis via isocyanate and secondary amine",
                "Urea synthesis via isocyanate and diazo",
                "Urea synthesis via isocyanate and sulfonamide",
                "urea",
            ]

            for reaction_type in urea_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Found urea disconnection strategy: {reaction_type}")
                    urea_disconnection_found = True
                    return

            # If no direct reaction match, check for the pattern manually
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for urea or thiourea in product
            product_has_urea = checker.check_fg("Urea", product)
            product_has_thiourea = checker.check_fg("Thiourea", product)

            if product_has_urea or product_has_thiourea:
                # Check for isocyanate/isothiocyanate and amine in reactants
                isocyanate_found = False
                amine_found = False

                for reactant in reactants:
                    if checker.check_fg("Isocyanate", reactant):
                        isocyanate_found = True
                    elif checker.check_fg("Isothiocyanate", reactant):
                        isocyanate_found = True

                    # Check for various types of amines
                    if (
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        or checker.check_fg("Aniline", reactant)
                    ):
                        amine_found = True

                if isocyanate_found and amine_found:
                    fg_type = "Urea" if product_has_urea else "Thiourea"
                    print(
                        f"Found {fg_type} disconnection strategy through pattern matching"
                    )
                    urea_disconnection_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return urea_disconnection_found
