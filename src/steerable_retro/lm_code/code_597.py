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
    This function detects if the synthesis involves an olefination reaction
    (conversion of aldehyde to terminal alkene).
    """
    olefination_present = False

    def dfs_traverse(node):
        nonlocal olefination_present

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # First check if this is a known olefination reaction type
            if (
                checker.check_reaction("Wittig reaction", rsmi)
                or checker.check_reaction("Julia Olefination", rsmi)
                or checker.check_reaction("Wittig with Phosphonium", rsmi)
            ):
                print(f"Known olefination reaction detected: {rsmi}")
                olefination_present = True
                return

            # If not a known reaction type, check for aldehyde to alkene conversion
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if any reactant has an aldehyde group
            reactant_has_aldehyde = False
            for reactant in reactants_part.split("."):
                if checker.check_fg("Aldehyde", reactant):
                    reactant_has_aldehyde = True
                    print(f"Reactant with aldehyde found: {reactant}")
                    break

            # Check if product has a terminal alkene
            product_has_terminal_alkene = False
            if checker.check_fg("Alkene", product_part) or checker.check_fg("Vinyl", product_part):
                product_has_terminal_alkene = True
                print(f"Product with terminal alkene found: {product_part}")

            # If both conditions are met, it's likely an olefination
            if reactant_has_aldehyde and product_has_terminal_alkene:
                print(f"Olefination from aldehyde detected at reaction: {rsmi}")
                olefination_present = True

        # Continue DFS traversal
        for child in node.get("children", []):
            if not olefination_present:  # Stop traversal if we already found an olefination
                dfs_traverse(child)

    dfs_traverse(route)
    print(f"Olefination from aldehyde present: {olefination_present}")
    return olefination_present
