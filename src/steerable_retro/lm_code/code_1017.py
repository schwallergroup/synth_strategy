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
    This function detects triazine ring formation from non-ring precursors.
    """
    triazine_formed = False

    def dfs_traverse(node):
        nonlocal triazine_formed

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product contains triazine but reactants don't
            product_mol = Chem.MolFromSmiles(product_smiles)
            triazine_pattern = Chem.MolFromSmarts("[n]1[c][n][c][n]1")

            if product_mol and triazine_pattern:
                product_has_triazine = product_mol.HasSubstructMatch(triazine_pattern)

                if product_has_triazine:
                    # Check if reactants don't have triazine
                    reactants_have_triazine = False
                    for r_smiles in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r_smiles)
                        if r_mol and r_mol.HasSubstructMatch(triazine_pattern):
                            reactants_have_triazine = True
                            break

                    if not reactants_have_triazine:
                        print(f"Triazine formation detected in reaction: {rsmi}")
                        triazine_formed = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return triazine_formed
