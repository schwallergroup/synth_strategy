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
    Detects if the synthesis route involves a late-stage ester hydrolysis
    (conversion of ester to carboxylic acid in the final step).
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # Check if this is a reaction node at a late stage (depth <= 1)
        if node["type"] == "reaction" and depth <= 1:
            try:
                # Extract reaction SMILES from metadata
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Checking late-stage reaction at depth {depth}: {rsmi}")

                    # Check if this is an ester hydrolysis reaction
                    is_ester_hydrolysis = checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                    )
                    print(f"Is ester hydrolysis reaction: {is_ester_hydrolysis}")

                    # If reaction checker failed, try checking functional groups manually
                    if not is_ester_hydrolysis:
                        # Check if any reactant contains an ester group
                        has_ester = any(
                            checker.check_fg("Ester", reactant) for reactant in reactants
                        )
                        print(f"Has ester in reactants: {has_ester}")

                        # Check if product contains a carboxylic acid group
                        has_acid = checker.check_fg("Carboxylic acid", product)
                        print(f"Has carboxylic acid in product: {has_acid}")

                        # Verify both conditions are met
                        if has_ester and has_acid:
                            print("Detected ester and carboxylic acid functional groups")
                            is_ester_hydrolysis = True

                    if is_ester_hydrolysis:
                        print(f"Detected late-stage ester hydrolysis at depth {depth}")
                        result = True
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    print("Starting traversal")
    dfs_traverse(route)
    print(f"Final result: {result}")
    return result
