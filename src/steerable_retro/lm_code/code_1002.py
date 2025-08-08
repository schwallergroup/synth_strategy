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
    Detects if the synthetic route builds complex nitrogen-rich heterocyclic structures.
    """
    nitrogen_heterocycle_count = 0
    morpholine_pattern = Chem.MolFromSmarts("[#7]1-[#6]-[#6]-[#8]-[#6]-[#6]-1")
    piperazine_pattern = Chem.MolFromSmarts("[#7]1-[#6]-[#6]-[#7]-[#6]-[#6]-1")
    piperidine_pattern = Chem.MolFromSmarts("[#7]1-[#6]-[#6]-[#6]-[#6]-[#6]-1")

    def dfs_traverse(node):
        nonlocal nitrogen_heterocycle_count

        if node["type"] == "mol" and node.get("smiles"):
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Count different nitrogen heterocycles
                    if mol.HasSubstructMatch(morpholine_pattern):
                        nitrogen_heterocycle_count += 1
                        print(f"Found morpholine structure in: {node['smiles']}")
                    if mol.HasSubstructMatch(piperazine_pattern):
                        nitrogen_heterocycle_count += 1
                        print(f"Found piperazine structure in: {node['smiles']}")
                    if mol.HasSubstructMatch(piperidine_pattern):
                        nitrogen_heterocycle_count += 1
                        print(f"Found piperidine structure in: {node['smiles']}")
            except:
                print(f"Error processing molecule SMILES: {node['smiles']}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total nitrogen heterocycles detected: {nitrogen_heterocycle_count}")
    return nitrogen_heterocycle_count >= 2
