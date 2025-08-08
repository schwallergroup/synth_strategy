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
    This function detects a strategy where a heterocyclic core is built and
    then modified through multiple steps, maintaining the core structure.
    """
    heterocycle_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Define patterns for various heterocycles
                thiophene_pattern = Chem.MolFromSmarts("c1cscc1")
                pyrimidine_pattern = Chem.MolFromSmarts("c1ncncn1")

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        if product_mol.HasSubstructMatch(
                            thiophene_pattern
                        ) or product_mol.HasSubstructMatch(pyrimidine_pattern):
                            heterocycle_depths.append(depth)
                            print(f"Heterocyclic core detected at depth {depth}")
                except:
                    print("Error processing heterocycle detection")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if heterocyclic core is maintained through multiple steps
    if len(heterocycle_depths) >= 3:  # Present in at least 3 steps
        print("Heterocyclic core building strategy detected")
        return True
    return False
