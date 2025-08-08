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


def main(route):
    """
    This function detects a linear synthesis strategy where each intermediate
    is used in only one subsequent reaction.
    """
    # Track molecules and their usage
    molecule_usage_count = {}
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Count usage of reactants
            for reactant in reactants_smiles:
                if reactant in molecule_usage_count:
                    molecule_usage_count[reactant] += 1
                    # If a molecule is used more than once, it's not a linear synthesis
                    if molecule_usage_count[reactant] > 1:
                        is_linear = False
                        print(f"Non-linear: molecule {reactant} used multiple times")
                else:
                    molecule_usage_count[reactant] = 1

            # Add product to tracking
            molecule_usage_count[product_smiles] = 0

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if synthesis is linear
    if is_linear:
        print("Detected linear synthesis strategy")

    return is_linear
