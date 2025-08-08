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
    This function detects isothiocyanate formation from amines.
    """
    isothiocyanate_formed = False

    def dfs_traverse(node):
        nonlocal isothiocyanate_formed

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            for r_smiles in reactants_smiles:
                r_mol = Chem.MolFromSmiles(r_smiles)
                if r_mol and amine_pattern and r_mol.HasSubstructMatch(amine_pattern):
                    # Check if product has isothiocyanate
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    isothiocyanate_pattern = Chem.MolFromSmarts("[N]=[C]=[S]")

                    if (
                        product_mol
                        and isothiocyanate_pattern
                        and product_mol.HasSubstructMatch(isothiocyanate_pattern)
                    ):
                        print(f"Isothiocyanate formation detected in reaction: {rsmi}")
                        isothiocyanate_formed = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return isothiocyanate_formed
