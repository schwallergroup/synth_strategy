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
    This function detects nitro reduction to amine in the synthetic route.
    """
    nitro_reduction_found = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro group in reactants
            nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")

            for r_smiles in reactants_smiles:
                r_mol = Chem.MolFromSmiles(r_smiles)
                if r_mol and nitro_pattern and r_mol.HasSubstructMatch(nitro_pattern):
                    # Check if product has amine at the same position
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    amine_pattern = Chem.MolFromSmarts("[NH2]")

                    if (
                        product_mol
                        and amine_pattern
                        and product_mol.HasSubstructMatch(amine_pattern)
                    ):
                        print(f"Nitro reduction detected in reaction: {rsmi}")
                        nitro_reduction_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return nitro_reduction_found
