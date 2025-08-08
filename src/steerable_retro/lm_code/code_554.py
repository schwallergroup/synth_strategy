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
    Detects a synthetic strategy involving multiple heterocyclic moieties
    (thiophene, benzimidazole, and indole).
    """
    has_thiophene = False
    has_benzimidazole = False
    has_indole = False

    def dfs_traverse(node):
        nonlocal has_thiophene, has_benzimidazole, has_indole

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for heterocycles
                    thiophene_pattern = Chem.MolFromSmarts("c1cscc1")
                    benzimidazole_pattern = Chem.MolFromSmarts("c1nc2ccccc2[nH]1")
                    indole_pattern = Chem.MolFromSmarts("c1[nH]c2ccccc2c1")

                    if mol.HasSubstructMatch(thiophene_pattern):
                        has_thiophene = True
                        print("Found thiophene moiety")

                    if mol.HasSubstructMatch(benzimidazole_pattern):
                        has_benzimidazole = True
                        print("Found benzimidazole moiety")

                    if mol.HasSubstructMatch(indole_pattern):
                        has_indole = True
                        print("Found indole moiety")
            except:
                print(f"Error processing molecule SMILES: {node['smiles']}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if at least two different heterocycles are found
    heterocycle_count = sum([has_thiophene, has_benzimidazole, has_indole])
    return heterocycle_count >= 2
