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
    This function detects a strategy where a diamine compound undergoes
    selective functionalization at one amine while the other remains intact
    or is protected/deprotected.
    """
    # Initialize tracking variables
    has_diamine = False
    has_selective_functionalization = False

    # SMARTS patterns
    diamine_pattern = Chem.MolFromSmarts("[NH2][c]1[cH][cH][cH][c]([NH2])[cH]1")
    sulfonamide_pattern = Chem.MolFromSmarts("[NH]S(=O)(=O)[c]")

    def dfs_traverse(node):
        nonlocal has_diamine, has_selective_functionalization

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(diamine_pattern):
                print("Detected diamine compound")
                has_diamine = True

        elif node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            product = Chem.MolFromSmiles(product_smiles)

            # Check for selective functionalization
            if product and product.HasSubstructMatch(sulfonamide_pattern):
                # Check if one amine is functionalized while another remains
                if product.GetSubstructMatches(Chem.MolFromSmarts("[NH2][c]")):
                    print("Detected selective functionalization of diamine")
                    has_selective_functionalization = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy was detected
    strategy_detected = has_diamine and has_selective_functionalization
    if strategy_detected:
        print("Detected diamine selective functionalization strategy")
    return strategy_detected
