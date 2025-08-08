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
    This function detects a synthesis strategy involving oxime chemistry.
    """
    oxime_chemistry_found = False

    def dfs_traverse(node):
        nonlocal oxime_chemistry_found

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            products_smiles = rsmi.split(">")[-1]

            # Check for oxime pattern
            oxime_pattern = Chem.MolFromSmarts("[#6]=[#7]-[#8]")

            reactant_mol = Chem.MolFromSmiles(reactants_smiles)
            product_mols = [Chem.MolFromSmiles(p) for p in products_smiles.split(".")]

            if oxime_pattern is not None:
                # Check if reactant or product has oxime
                if reactant_mol is not None and reactant_mol.HasSubstructMatch(oxime_pattern):
                    print("Found oxime in reactants")
                    oxime_chemistry_found = True

                for p_mol in product_mols:
                    if p_mol is not None and p_mol.HasSubstructMatch(oxime_pattern):
                        print("Found oxime in products")
                        oxime_chemistry_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return oxime_chemistry_found
