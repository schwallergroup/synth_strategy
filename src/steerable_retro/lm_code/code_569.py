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
    Detects a synthetic strategy involving multiple halogenated aromatics.
    Looks for the presence of at least two different halogens (F, Cl, Br, I) on aromatic rings.
    """
    # Initialize counters for different halogens
    has_fluoro = False
    has_chloro = False
    has_bromo = False
    has_iodo = False

    # SMARTS patterns for halogenated aromatics
    fluoro_pattern = Chem.MolFromSmarts("[c][F]")
    chloro_pattern = Chem.MolFromSmarts("[c][Cl]")
    bromo_pattern = Chem.MolFromSmarts("[c][Br]")
    iodo_pattern = Chem.MolFromSmarts("[c][I]")

    def dfs_traverse(node):
        nonlocal has_fluoro, has_chloro, has_bromo, has_iodo

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if mol.HasSubstructMatch(fluoro_pattern):
                    has_fluoro = True
                    print("Found fluorinated aromatic")
                if mol.HasSubstructMatch(chloro_pattern):
                    has_chloro = True
                    print("Found chlorinated aromatic")
                if mol.HasSubstructMatch(bromo_pattern):
                    has_bromo = True
                    print("Found brominated aromatic")
                if mol.HasSubstructMatch(iodo_pattern):
                    has_iodo = True
                    print("Found iodinated aromatic")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Count the number of different halogens present
    halogen_count = sum([has_fluoro, has_chloro, has_bromo, has_iodo])

    # Strategy is present if at least two different halogens are found
    strategy_present = halogen_count >= 2

    if strategy_present:
        print(
            f"Detected multi-halogenated aromatics strategy with {halogen_count} different halogens"
        )
    else:
        print(f"Multi-halogenated strategy not detected. Found only {halogen_count} halogen types")

    return strategy_present
