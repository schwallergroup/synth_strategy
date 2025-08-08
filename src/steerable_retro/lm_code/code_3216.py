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
    This function detects if the synthesis uses a linear strategy with
    multiple ether formations, particularly focusing on C-O bond formations.
    """
    # Track ether formations
    ether_formations = 0
    is_linear = True  # Assume linear until proven otherwise

    def dfs_traverse(node):
        nonlocal ether_formations, is_linear

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]
            reactants = reactants_part.split(".")
            product = product_part

            # Check if this is a convergent step (more than 2 reactants)
            if len(reactants) > 2:
                is_linear = False
                print("Found convergent step with more than 2 reactants")

            # Check for ether formation
            ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]")

            product_mol = Chem.MolFromSmiles(product) if product else None
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

            if product_mol:
                product_ethers = len(product_mol.GetSubstructMatches(ether_pattern))
                reactant_ethers = sum(
                    len(r.GetSubstructMatches(ether_pattern)) for r in reactant_mols if r
                )

                if product_ethers > reactant_ethers:
                    ether_formations += 1
                    print(f"Found ether formation, total: {ether_formations}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if it's a linear synthesis with at least 2 ether formations
    return is_linear and ether_formations >= 2
