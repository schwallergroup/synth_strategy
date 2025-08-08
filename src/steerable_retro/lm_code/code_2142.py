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
    Detects if the synthesis uses a late-stage aldol condensation (depth 0-1)
    for C-C bond formation between a carbonyl and an activated methylene.
    """
    aldol_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal aldol_detected

        if node["type"] == "reaction" and depth <= 1:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for aldehyde pattern in reactants
            aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
            activated_methylene_pattern = Chem.MolFromSmarts("[CX4H2]([C]=[O,N])")

            # Check for C=C-C=O pattern in product (aldol condensation result)
            aldol_product_pattern = Chem.MolFromSmarts("[C]=[C]-[C]=[O,N]")

            has_aldehyde = False
            has_activated_methylene = False

            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(aldehyde_pattern):
                    has_aldehyde = True
                if mol and mol.HasSubstructMatch(activated_methylene_pattern):
                    has_activated_methylene = True

            product_mol = Chem.MolFromSmiles(product_smiles)
            has_aldol_product = product_mol and product_mol.HasSubstructMatch(aldol_product_pattern)

            if has_aldehyde and has_activated_methylene and has_aldol_product:
                print(f"Late-stage aldol condensation detected at depth {depth}")
                aldol_detected = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return aldol_detected
