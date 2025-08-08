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
    This function detects a Boc protection-deprotection sequence in the synthesis.
    """
    boc_protected_found = False
    deprotection_found = False

    def dfs_traverse(node):
        nonlocal boc_protected_found, deprotection_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc group in molecules
            boc_pattern = Chem.MolFromSmarts("C(=O)OC(C)(C)C")

            # Check for deprotection (Boc in reactants but not in product)
            reactants_have_boc = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(boc_pattern) for r in reactants if r
            )
            product_has_boc = (
                Chem.MolFromSmiles(product).HasSubstructMatch(boc_pattern) if product else False
            )

            if reactants_have_boc and not product_has_boc:
                deprotection_found = True
                print("Boc deprotection detected")

            # Check for protection (Boc in product)
            if product_has_boc:
                boc_protected_found = True
                print("Boc protected intermediate detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both protection and deprotection were found
    result = boc_protected_found and deprotection_found
    print(f"Boc protection-deprotection sequence detected: {result}")
    return result
