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
    Detects if the synthesis follows a linear strategy where each reaction
    has only one product that serves as a reactant in the next step.
    """
    is_linear = True
    reaction_count = 0

    # Track the synthesis path
    reaction_nodes = []

    def dfs_traverse(node):
        nonlocal is_linear, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1
            reaction_nodes.append(node)

            # Check if this reaction has multiple products
            try:
                rsmi = node["metadata"]["rsmi"]
                products_part = rsmi.split(">")[-1]
                products = [p for p in products_part.split(".") if p.strip()]

                if len(products) > 1:
                    # If a reaction produces multiple significant products, it's not strictly linear
                    print(f"Found non-linear step with {len(products)} products")
                    is_linear = False
            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")
                is_linear = False

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # A route with only one reaction is trivially linear
    if reaction_count <= 1:
        return True

    # Check connectivity between reactions
    # In a linear synthesis, each reaction should have exactly one child reaction
    # except for the last reaction in the sequence
    for node in route.get("children", []):
        if node["type"] == "reaction":
            reaction_children = [
                child for child in node.get("children", []) if child["type"] == "reaction"
            ]

            # If this reaction has more than one child reaction, it's not linear
            if len(reaction_children) > 1:
                print(f"Found reaction with multiple child reactions: {len(reaction_children)}")
                is_linear = False
                break

    return is_linear
