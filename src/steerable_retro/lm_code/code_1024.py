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
    Detects the use of Boc protection/deprotection in the synthetic route.
    """
    boc_protection_found = False

    def dfs_traverse(node):
        nonlocal boc_protection_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Check for Boc group in product
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol:
                boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)[N]")
                if product_mol.HasSubstructMatch(boc_pattern):
                    # Check if Boc was introduced in this step
                    boc_in_reactants = False
                    for r in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol and r_mol.HasSubstructMatch(boc_pattern):
                            boc_in_reactants = True
                            break

                    if not boc_in_reactants:
                        boc_protection_found = True
                        print(f"Found Boc protection: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return boc_protection_found
