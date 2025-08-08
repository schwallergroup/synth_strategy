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
    Detects the use of silyl protecting groups in the synthesis.
    """
    has_silyl_protection = False

    def dfs_traverse(node):
        nonlocal has_silyl_protection

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for silyl group pattern
            silyl_pattern = Chem.MolFromSmarts("[#6]-[#14](-[#6])(-[#6])-[#6]")

            # Check reactants or products for silyl groups
            for smiles in reactants_smiles + [product_smiles]:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol and mol.HasSubstructMatch(silyl_pattern):
                        has_silyl_protection = True
                        print("Detected silyl protecting group")
                        return
                except:
                    continue

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return has_silyl_protection
