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
    Detects if the synthesis uses Boc deprotection as part of a protection/deprotection strategy.
    """
    boc_deprotection_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_deprotection_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for Boc group in reactant
            boc_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[OX2][CX4]([CX4])([CX4])[CX4]")

            # Check for free amine in product
            amine_pattern = Chem.MolFromSmarts("[NX3;H1,H2]")

            has_boc = False
            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(boc_pattern):
                    has_boc = True
                    break

            product_mol = Chem.MolFromSmiles(product_smiles)
            has_amine = product_mol and product_mol.HasSubstructMatch(amine_pattern)

            # Check if Boc is removed (present in reactant but not in product)
            if has_boc and has_amine and not product_mol.HasSubstructMatch(boc_pattern):
                print(f"Boc deprotection detected at depth {depth}")
                boc_deprotection_detected = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return boc_deprotection_detected
