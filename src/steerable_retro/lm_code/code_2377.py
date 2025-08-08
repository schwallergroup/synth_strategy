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
    Detects the use of halogenated aromatic compounds throughout the synthesis.
    """
    halogen_count = 0

    def dfs_traverse(node):
        nonlocal halogen_count

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # Check for halogenated aromatics
                fluoro_aromatic = Chem.MolFromSmarts("[F][c]")
                chloro_aromatic = Chem.MolFromSmarts("[Cl][c]")
                bromo_aromatic = Chem.MolFromSmarts("[Br][c]")

                if mol.HasSubstructMatch(fluoro_aromatic):
                    halogen_count += len(mol.GetSubstructMatches(fluoro_aromatic))
                if mol.HasSubstructMatch(chloro_aromatic):
                    halogen_count += len(mol.GetSubstructMatches(chloro_aromatic))
                if mol.HasSubstructMatch(bromo_aromatic):
                    halogen_count += len(mol.GetSubstructMatches(bromo_aromatic))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Found {halogen_count} halogenated aromatic positions")
    return halogen_count >= 3  # At least 3 halogen atoms on aromatic rings
