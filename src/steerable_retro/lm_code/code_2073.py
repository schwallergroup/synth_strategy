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
    Detects a synthetic route that includes tetrahydropyran (THP) deprotection of a nitrogen.
    """
    has_thp_deprotection = False

    def dfs_traverse(node):
        nonlocal has_thp_deprotection

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for THP deprotection
            thp_pattern = Chem.MolFromSmarts("[#7]C1CCCCO1")
            nh_pattern = Chem.MolFromSmarts("[#7;H]")

            if any(
                mol.HasSubstructMatch(thp_pattern) for mol in reactants
            ) and product.HasSubstructMatch(nh_pattern):
                print("Found THP deprotection")
                has_thp_deprotection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"THP deprotection strategy detected: {has_thp_deprotection}")
    return has_thp_deprotection
