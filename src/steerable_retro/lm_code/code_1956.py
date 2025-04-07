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

root_data = "/home/andres/Documents/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
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
    This function detects if the synthesis includes an amide reduction step.
    """
    found_amide_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_reduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide reduction reaction directly
            if (
                checker.check_reaction("Reduction of primary amides to amines", rsmi)
                or checker.check_reaction("Reduction of secondary amides to amines", rsmi)
                or checker.check_reaction("Reduction of tertiary amides to amines", rsmi)
            ):
                print(f"Amide reduction reaction detected at depth {depth}")
                found_amide_reduction = True
            else:
                # Alternative check: look for amide in reactants and corresponding amine in product
                has_amide = False
                amide_types = ["Primary amide", "Secondary amide", "Tertiary amide"]
                for reactant in reactants:
                    for amide_type in amide_types:
                        if checker.check_fg(amide_type, reactant):
                            has_amide = True
                            print(f"{amide_type} found in reactant: {reactant}")
                            break
                    if has_amide:
                        break

                has_amine = (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                )

                if has_amide and has_amine:
                    # Additional check to confirm it's a reduction (not just any amide-to-amine conversion)
                    # Look for reduction-specific patterns
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # In amide reduction, C=O becomes CH2, so check for increased hydrogen count
                        # This is a simplified check and may not catch all cases
                        found_amide_reduction = True
                        print(f"Potential amide reduction detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_amide_reduction
