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
    This function detects if the synthesis uses Boc protection of an amine
    as part of the synthetic strategy.
    """
    boc_protection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for primary amine in reactants
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            # Check for Boc-protected amine in product
            boc_pattern = Chem.MolFromSmarts("[NH]C(=O)OC(C)(C)C")

            # Check if transformation is amine to Boc-protected amine
            amine_present = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(amine_pattern):
                    amine_present = True
                    break

            product_mol = Chem.MolFromSmiles(product)
            boc_present = product_mol and product_mol.HasSubstructMatch(boc_pattern)

            if amine_present and boc_present:
                print(f"Boc protection detected at depth {depth}")
                boc_protection_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return boc_protection_found
