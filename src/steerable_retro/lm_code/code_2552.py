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
    This function detects a synthesis strategy involving alcohol oxidation/esterification.
    """
    alcohol_oxidation_found = False

    def dfs_traverse(node):
        nonlocal alcohol_oxidation_found

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            products_smiles = rsmi.split(">")[-1]

            # Check for alcohol to ester transformation
            alcohol_pattern = Chem.MolFromSmarts("[#6][#8H]")
            ester_pattern = Chem.MolFromSmarts("[#6](=[#8])[#8][#6]")

            reactant_mol = Chem.MolFromSmiles(reactants_smiles)
            product_mol = Chem.MolFromSmiles(products_smiles)

            if (
                reactant_mol is not None
                and product_mol is not None
                and alcohol_pattern is not None
                and ester_pattern is not None
            ):
                if reactant_mol.HasSubstructMatch(
                    alcohol_pattern
                ) and product_mol.HasSubstructMatch(ester_pattern):
                    print("Found alcohol oxidation/esterification")
                    alcohol_oxidation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return alcohol_oxidation_found
