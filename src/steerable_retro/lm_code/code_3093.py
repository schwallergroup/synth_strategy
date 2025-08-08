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
    This function detects a strategy where a nitrile intermediate is used
    and later incorporated into heterocycle formation.
    """
    nitrile_found = False
    heterocycle_after_nitrile = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_found, heterocycle_after_nitrile

        if node["type"] == "mol":
            # Check for nitrile group
            if node.get("smiles"):
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    nitrile_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
                    if mol.HasSubstructMatch(nitrile_pattern):
                        print(f"Found nitrile at depth {depth}")
                        nitrile_found = True

        elif node["type"] == "reaction" and nitrile_found:
            # Check if product contains heterocycle after nitrile was found
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]
                mol = Chem.MolFromSmiles(product)
                if mol:
                    isoxazole_pattern = Chem.MolFromSmarts("c1onc(c1)")
                    pyrazole_pattern = Chem.MolFromSmarts("c1cnn(c1)")

                    if mol.HasSubstructMatch(isoxazole_pattern) or mol.HasSubstructMatch(
                        pyrazole_pattern
                    ):
                        print(f"Found heterocycle formation after nitrile at depth {depth}")
                        heterocycle_after_nitrile = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return nitrile_found and heterocycle_after_nitrile
