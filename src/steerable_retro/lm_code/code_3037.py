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
    This function detects a linear synthesis strategy where each reaction
    has only one product that becomes a reactant in the next step.
    """
    is_linear = True
    reaction_count = 0
    intermediate_molecules = set()  # Track intermediate molecules

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1

            # Check if this reaction has exactly one product
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                products = rsmi.split(">")[-1].split(".")

                # If more than one product, it's not a strictly linear synthesis
                if len(products) > 1:
                    is_linear = False
                    print(f"Non-linear reaction detected at depth {depth}: multiple products")

            # In a retrosynthetic tree, a reaction node should have 1 or 2 children (reactants)
            # If it has 0 or >2 children, it's not a linear synthesis
            if len(node.get("children", [])) == 0 or len(node.get("children", [])) > 2:
                is_linear = False
                print(
                    f"Non-linear reaction detected at depth {depth}: {len(node.get('children', []))} reactants"
                )

            # For each child molecule, mark it as an intermediate
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    intermediate_molecules.add(child["smiles"])

        elif node["type"] == "mol" and not node.get("in_stock", False):
            # For intermediate molecules, check if they participate in exactly one reaction
            # In a linear synthesis, an intermediate molecule should have at most one reaction child
            reaction_children = [
                child for child in node.get("children", []) if child["type"] == "reaction"
            ]
            if len(reaction_children) > 1:
                is_linear = False
                print(
                    f"Non-linear synthesis detected at depth {depth}: molecule participates in multiple reactions"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # A linear synthesis should have at least 2 reactions
    if reaction_count < 2:
        is_linear = False
        print(f"Not enough reactions for linear synthesis: {reaction_count}")

    if is_linear:
        print(f"Detected linear synthesis strategy with {reaction_count} reactions")
    else:
        print("Linear synthesis strategy not detected")

    return is_linear
