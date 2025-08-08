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
    This function detects the conversion of a cyano group to an amide group.
    """
    cn_to_amide_conversion_detected = False

    def dfs_traverse(node):
        nonlocal cn_to_amide_conversion_detected

        if node["type"] == "reaction" and not cn_to_amide_conversion_detected:
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Examining reaction: {rsmi}")

                # First check if this is a known nitrile to amide reaction
                if checker.check_reaction("Nitrile to amide", rsmi) or checker.check_reaction(
                    "Nitrile and hydrogen peroxide to amide", rsmi
                ):
                    print(f"Known nitrile to amide reaction detected: {rsmi}")
                    cn_to_amide_conversion_detected = True
                    return

                # If not a known reaction type, check for functional group conversion
                for reactant_smiles in reactants_smiles:
                    # Check if reactant contains nitrile group
                    if checker.check_fg("Nitrile", reactant_smiles):
                        print(f"Reactant contains nitrile: {reactant_smiles}")

                        # Check if product contains amide group
                        if (
                            checker.check_fg("Primary amide", product_smiles)
                            or checker.check_fg("Secondary amide", product_smiles)
                            or checker.check_fg("Tertiary amide", product_smiles)
                        ):

                            print(f"Product contains amide: {product_smiles}")
                            print(f"Cyano to amide conversion detected in reaction: {rsmi}")
                            cn_to_amide_conversion_detected = True
                            return
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return cn_to_amide_conversion_detected
