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
    Detects if the synthetic route involves heterocyclic compounds containing methoxy groups.
    """
    has_methoxy_heterocycle = False

    def dfs_traverse(node):
        nonlocal has_methoxy_heterocycle

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Pattern for methoxy group
                methoxy_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]")
                # Pattern for heterocycle (simplified)
                heterocycle_pattern = Chem.MolFromSmarts(
                    "[#6]1[#6,#7,#8,#16][#6,#7,#8,#16][#6,#7,#8,#16][#6,#7,#8,#16][#6,#7,#8,#16]1"
                )

                if mol.HasSubstructMatch(methoxy_pattern) and mol.HasSubstructMatch(
                    heterocycle_pattern
                ):
                    has_methoxy_heterocycle = True
                    print(f"Found methoxy-containing heterocycle in molecule: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_methoxy_heterocycle
