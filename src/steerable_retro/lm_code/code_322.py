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
    This function detects if a sulfonamide and chloroaryl moiety are maintained
    throughout the synthesis route.
    """
    # Track if these functional groups are present at each step
    all_products_have_sulfonamide = True
    all_products_have_chloroaryl = True

    def dfs_traverse(node, depth=0):
        nonlocal all_products_have_sulfonamide, all_products_have_chloroaryl

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]

            try:
                product = Chem.MolFromSmiles(product_smiles)

                # Check for sulfonamide group
                sulfonamide_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])[#7]")
                if not (product and product.HasSubstructMatch(sulfonamide_pattern)):
                    all_products_have_sulfonamide = False
                    print(f"Product at depth {depth} lacks sulfonamide group")

                # Check for chloroaryl moiety
                chloroaryl_pattern = Chem.MolFromSmarts("[c]-[Cl]")
                if not (product and product.HasSubstructMatch(chloroaryl_pattern)):
                    all_products_have_chloroaryl = False
                    print(f"Product at depth {depth} lacks chloroaryl moiety")

            except:
                # If we can't parse the SMILES, assume the pattern isn't maintained
                all_products_have_sulfonamide = False
                all_products_have_chloroaryl = False

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return all_products_have_sulfonamide and all_products_have_chloroaryl
