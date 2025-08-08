#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    Detects conversion of an ester to an alcohol
    """
    conversion_found = False

    def dfs_traverse(node, depth=0):
        nonlocal conversion_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Depth {depth}, Examining reaction: {rsmi}")

            # Check if any reactant contains an ester group
            reactant_has_ester = any(checker.check_fg("Ester", r) for r in reactants)

            # Check if product contains an alcohol group
            product_has_alcohol = (
                checker.check_fg("Primary alcohol", product)
                or checker.check_fg("Secondary alcohol", product)
                or checker.check_fg("Tertiary alcohol", product)
            )

            # Check if this is a reduction reaction
            is_reduction_reaction = checker.check_reaction(
                "Reduction of ester to primary alcohol", rsmi
            )

            print(f"  Reactant has ester: {reactant_has_ester}")
            print(f"  Product has alcohol: {product_has_alcohol}")
            print(f"  Is reduction reaction: {is_reduction_reaction}")

            # If direct reduction reaction is detected
            if is_reduction_reaction:
                conversion_found = True
                print("  Found ester to alcohol conversion via reduction reaction")
            # Fallback check for ester to alcohol conversion
            elif reactant_has_ester and product_has_alcohol:
                # This is a more general check that might catch other reaction types
                # that convert esters to alcohols
                conversion_found = True
                print("  Found ester to alcohol conversion")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return conversion_found
