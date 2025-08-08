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
    Detects if the synthesis involves ester aminolysis to form an amide.
    """
    ester_aminolysis_detected = False

    def dfs_traverse(node):
        nonlocal ester_aminolysis_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # First, check if this is an ester aminolysis reaction directly
                if checker.check_reaction("Aminolysis of esters", rsmi):
                    print(f"Ester aminolysis detected directly: {rsmi}")
                    ester_aminolysis_detected = True
                    return

                # If not directly identified, check for the pattern manually
                if len(reactants) >= 2:
                    # Check for ester in reactants
                    ester_reactant = None
                    amine_reactant = None

                    for reactant in reactants:
                        if checker.check_fg("Ester", reactant):
                            ester_reactant = reactant
                        if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                            "Secondary amine", reactant
                        ):
                            amine_reactant = reactant

                    # Check if product has an amide
                    if (
                        ester_reactant
                        and amine_reactant
                        and checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                    ):
                        # Verify this is not just a coincidence by checking the reaction pattern
                        # This is a fallback in case the direct reaction check failed
                        print(f"Ester aminolysis pattern detected: {rsmi}")
                        print(f"Ester reactant: {ester_reactant}")
                        print(f"Amine reactant: {amine_reactant}")
                        print(f"Amide product: {product}")
                        ester_aminolysis_detected = True
            except Exception as e:
                print(f"Error in processing reaction SMILES: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return ester_aminolysis_detected
