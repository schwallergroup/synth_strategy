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
    This function detects a synthetic strategy involving halogenated aromatic rings
    (chloro- and fluoro-substituted) in the final product.
    """
    has_halogenated_aromatics = False

    def dfs_traverse(node, depth=0):
        nonlocal has_halogenated_aromatics

        if node["type"] == "mol" and depth == 0:  # Final product
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # Check for chloro-aromatic pattern
                cl_aromatic_pattern = Chem.MolFromSmarts("c[Cl]")
                # Check for fluoro-aromatic pattern
                f_aromatic_pattern = Chem.MolFromSmarts("c[F]")

                has_cl = mol.HasSubstructMatch(cl_aromatic_pattern)
                has_f = mol.HasSubstructMatch(f_aromatic_pattern)

                if has_cl and has_f:
                    print("Found both chloro- and fluoro-substituted aromatic rings")
                    has_halogenated_aromatics = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_halogenated_aromatics
