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
    Detects the use of acetal protection for an aldehyde.
    Looks for aldehyde → acetal transformation.
    """
    acetal_protection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal acetal_protection_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aldehyde in reactants
            aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
            aldehyde_present = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(aldehyde_pattern):
                    aldehyde_present = True

            # Check for acetal in product
            acetal_pattern = Chem.MolFromSmarts("[#6]1[#8][#6][#6][#8]1")
            product_mol = Chem.MolFromSmiles(product)

            if aldehyde_present and product_mol and product_mol.HasSubstructMatch(acetal_pattern):
                acetal_protection_found = True
                print(f"Found acetal protection at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return acetal_protection_found
