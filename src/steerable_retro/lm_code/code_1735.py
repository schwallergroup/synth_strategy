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
    This function detects if the synthesis involves multiple types of heterocycles.
    """
    heterocycle_types = set()

    def dfs_traverse(node):
        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for different heterocycle patterns
                indole_pattern = Chem.MolFromSmarts("[nH]1c2ccccc2cc1")
                pyridine_pattern = Chem.MolFromSmarts("n1ccccc1")
                phthalimide_pattern = Chem.MolFromSmarts("N1C(=O)c2ccccc2C1=O")

                if mol.HasSubstructMatch(indole_pattern):
                    heterocycle_types.add("indole")
                if mol.HasSubstructMatch(pyridine_pattern):
                    heterocycle_types.add("pyridine")
                if mol.HasSubstructMatch(phthalimide_pattern):
                    heterocycle_types.add("phthalimide")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If we have at least 3 different heterocycle types, it's a heterocycle-rich synthesis
    if len(heterocycle_types) >= 3:
        print(f"Heterocycle-rich synthesis detected with types: {', '.join(heterocycle_types)}")
        return True

    return False
