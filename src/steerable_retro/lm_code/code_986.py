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
    This function detects if the synthesis involves a beta-lactam core structure
    that is preserved throughout the synthesis.
    """
    # Beta-lactam is a 4-membered ring with N and C=O
    # Define the beta-lactam pattern correctly
    beta_lactam_pattern = Chem.MolFromSmarts("[#6]1[#6](=[O])[#7][#6]1")

    # Track if we've found a beta-lactam that's preserved
    beta_lactam_preserved = [False]

    def check_beta_lactam(smiles):
        """Helper function to check if a molecule contains a beta-lactam"""
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.HasSubstructMatch(beta_lactam_pattern):
            return True
        return False

    def dfs_traverse(node, depth=0):
        # Check molecule nodes for beta-lactam
        if node["type"] == "mol" and "smiles" in node:
            if check_beta_lactam(node["smiles"]):
                print(f"Beta-lactam found in molecule at depth {depth}: {node['smiles']}")

                # If this is the final product (depth 0), mark as preserved initially
                if depth == 0:
                    beta_lactam_preserved[0] = True

        # Check reaction nodes to see if beta-lactam is preserved
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_has_beta_lactam = check_beta_lactam(product)
                reactant_has_beta_lactam = any(check_beta_lactam(r) for r in reactants)

                # If the product has a beta-lactam but no reactant does, it's formed in this step
                # If a reactant has a beta-lactam but the product doesn't, it's destroyed
                if product_has_beta_lactam and not reactant_has_beta_lactam:
                    print(f"Beta-lactam formed in reaction at depth {depth}: {rsmi}")
                    beta_lactam_preserved[0] = False
                elif not product_has_beta_lactam and reactant_has_beta_lactam:
                    print(f"Beta-lactam destroyed in reaction at depth {depth}: {rsmi}")
                    beta_lactam_preserved[0] = False
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Continue traversing children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return beta_lactam_preserved[0]
