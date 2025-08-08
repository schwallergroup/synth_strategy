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
    Detects a synthetic strategy involving piperazine as a linker to heterocycles.
    """
    piperazine_present = False
    piperazine_linked_to_heterocycle = False

    def dfs_traverse(node):
        nonlocal piperazine_present, piperazine_linked_to_heterocycle

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for piperazine
                piperazine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6][#6]1")
                if mol.HasSubstructMatch(piperazine_pattern):
                    piperazine_present = True

                # Check for piperazine linked to heterocycle
                piperazine_linked_pattern = Chem.MolFromSmarts(
                    "[#7]1[#6][#6][#7]([#6]2[#7][#6][#6][#6][#6]2)[#6][#6]1"
                )
                if mol.HasSubstructMatch(piperazine_linked_pattern):
                    piperazine_linked_to_heterocycle = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Piperazine present: {piperazine_present}")
    print(f"Piperazine linked to heterocycle: {piperazine_linked_to_heterocycle}")

    return piperazine_present and piperazine_linked_to_heterocycle
