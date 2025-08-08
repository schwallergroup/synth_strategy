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
    This function detects if the synthesis route involves the formation of a phthalimide-like structure.
    """
    has_phthalimide_formation = False

    def dfs_traverse(node):
        nonlocal has_phthalimide_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]

            # Check if product contains phthalimide-like structure
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol:
                phthalimide_pattern = Chem.MolFromSmarts("[N]1C(=O)c2ccccc2C1=O")
                if product_mol.HasSubstructMatch(phthalimide_pattern):
                    has_phthalimide_formation = True
                    print(f"Detected phthalimide formation in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_phthalimide_formation
