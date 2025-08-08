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
    Detects if the synthesis involves creation of two or more stereocenters in a single step
    """
    stereoselective_step_found = False

    def dfs_traverse(node):
        nonlocal stereoselective_step_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Count stereocenters in reactants
            reactant_stereocenters = 0
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Count chiral centers using '@' symbol in SMILES as a simple heuristic
                        reactant_stereocenters += reactant.count("@")
                except:
                    continue

            # Count stereocenters in product
            product_stereocenters = 0
            try:
                product_stereocenters = product.count("@")
            except:
                pass

            # Check if at least 2 new stereocenters were created
            if product_stereocenters >= reactant_stereocenters + 2:
                print(
                    f"Stereoselective transformation detected: {reactant_stereocenters} -> {product_stereocenters} stereocenters"
                )
                stereoselective_step_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return stereoselective_step_found
