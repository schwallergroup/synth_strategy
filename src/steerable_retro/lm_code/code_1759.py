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
    Detects if the synthesis involves multiple halogen atoms (particularly fluorine)
    on aromatic rings in the final product.
    """
    # Initialize tracking variable
    multiple_halogens = False

    def dfs_traverse(node, depth=0):
        nonlocal multiple_halogens

        if node["type"] == "mol" and depth == 0:  # Final product
            # Check for multiple fluorine atoms on aromatic rings
            fluoro_aromatic_pattern = Chem.MolFromSmarts("[c]-[F]")

            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Count matches
                    matches = mol.GetSubstructMatches(fluoro_aromatic_pattern)
                    if len(matches) >= 2:
                        multiple_halogens = True
                        print(
                            f"Detected {len(matches)} fluorine atoms on aromatic rings in final product"
                        )
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    print(f"Halogenated aromatic strategy detected: {multiple_halogens}")
    return multiple_halogens
