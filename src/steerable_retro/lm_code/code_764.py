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
    Detects if the route has a late-stage ether formation (depth 0-1).
    """
    late_stage_ether_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_ether_formation

        if node["type"] == "reaction" and depth <= 1:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check for ether formation
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
            product = Chem.MolFromSmiles(product_smiles)

            # Count ethers in reactants and product
            ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]")

            reactant_ether_count = sum(
                [len(r.GetSubstructMatches(ether_pattern)) for r in reactants if r]
            )
            product_ether_count = len(product.GetSubstructMatches(ether_pattern)) if product else 0

            if product_ether_count > reactant_ether_count:
                print(f"Late-stage ether formation detected at depth {depth}")
                late_stage_ether_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return late_stage_ether_formation
