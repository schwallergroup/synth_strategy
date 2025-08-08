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
    This function detects a synthetic strategy involving multiple protection steps
    (both amine and aldehyde protection).
    """
    # Track if we found the key features
    found_amine_protection = False
    found_aldehyde_protection = False

    # SMARTS patterns
    boc_pattern = Chem.MolFromSmarts("[#6]C([#6])([#6])[#8]C(=[#8])[#7]")
    acetal_pattern = Chem.MolFromSmarts("[#6]1[#8][#6][#6][#8]1")

    def dfs_traverse(node, depth=0):
        nonlocal found_amine_protection, found_aldehyde_protection

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")

            if not rsmi:
                return

            reactants, products = rsmi.split(">")[0], rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(products)

                # Check for Boc protection
                if product_mol and product_mol.HasSubstructMatch(boc_pattern):
                    print(f"Found Boc protection at depth {depth}")
                    found_amine_protection = True

                # Check for acetal protection
                if product_mol and product_mol.HasSubstructMatch(acetal_pattern):
                    print(f"Found acetal protection at depth {depth}")
                    found_aldehyde_protection = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if both protection types were found
    return found_amine_protection and found_aldehyde_protection
