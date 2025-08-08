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
    Detects if the route has an early biaryl coupling (depth > 5).
    """
    early_biaryl_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal early_biaryl_coupling

        if node["type"] == "reaction" and depth > 5:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check for biaryl formation
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
            product = Chem.MolFromSmiles(product_smiles)

            # Biaryl pattern
            biaryl_pattern = Chem.MolFromSmarts("c-c")

            reactant_biaryl_count = sum(
                [len(r.GetSubstructMatches(biaryl_pattern)) for r in reactants if r]
            )
            product_biaryl_count = (
                len(product.GetSubstructMatches(biaryl_pattern)) if product else 0
            )

            if product_biaryl_count > reactant_biaryl_count:
                print(f"Early biaryl coupling detected at depth {depth}")
                early_biaryl_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return early_biaryl_coupling
