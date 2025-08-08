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
    This function detects a synthetic strategy involving selective functionalization
    of a symmetrical diol to create an asymmetric molecule.
    """
    diol_used = False
    asymmetric_product = False

    # SMARTS for diol
    diol_pattern = Chem.MolFromSmarts("[OH][CH2][CH2][CH2][CH2][CH2][OH]")

    def dfs_traverse(node):
        nonlocal diol_used, asymmetric_product

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                products_smiles = rsmi.split(">")[-1]

                try:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                    product = Chem.MolFromSmiles(products_smiles)

                    # Check if a diol is used as reactant
                    if any(r and r.HasSubstructMatch(diol_pattern) for r in reactants if r):
                        diol_used = True
                        print("Diol detected as reactant")

                    # Check if product has one OH and one protected or functionalized OH
                    if product:
                        oh_pattern = Chem.MolFromSmarts("[OH]")
                        oh_count = len(product.GetSubstructMatches(oh_pattern))

                        # If product has exactly one OH and originally came from a diol
                        if oh_count == 1 and diol_used:
                            asymmetric_product = True
                            print("Asymmetric functionalization of diol detected")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return diol_used and asymmetric_product
