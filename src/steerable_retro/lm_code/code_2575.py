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
    Detects if a beta-lactam ring is preserved throughout the synthesis
    and present in the final product.
    """
    # Track beta-lactam presence
    has_beta_lactam_in_final = False
    has_beta_lactam_in_intermediates = False

    # SMARTS pattern for beta-lactam
    beta_lactam_pattern = "[#6]1[#6][#7][#6](=[#8])1"

    def dfs_traverse(node, depth=0):
        nonlocal has_beta_lactam_in_final, has_beta_lactam_in_intermediates

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(Chem.MolFromSmarts(beta_lactam_pattern)):
                if depth == 0:  # Final product
                    has_beta_lactam_in_final = True
                    print("Detected beta-lactam in final product")
                else:  # Intermediate
                    has_beta_lactam_in_intermediates = True
                    print(f"Detected beta-lactam in intermediate at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    preserves_beta_lactam = has_beta_lactam_in_final and has_beta_lactam_in_intermediates

    if preserves_beta_lactam:
        print("Detected beta-lactam preservation strategy")

    return preserves_beta_lactam
