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
    This function detects the use of fluorinated building blocks in the synthesis,
    specifically trifluoromethyl sulfide and trifluorovinyl groups.
    """
    # Initialize flags
    has_trifluoromethyl_sulfide = False
    has_trifluorovinyl = False

    def dfs_traverse(node):
        nonlocal has_trifluoromethyl_sulfide, has_trifluorovinyl

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for trifluoromethyl sulfide group
                trifluoromethyl_sulfide_pattern = Chem.MolFromSmarts("[#16][C]([F])([F])[F]")
                if mol.HasSubstructMatch(trifluoromethyl_sulfide_pattern):
                    has_trifluoromethyl_sulfide = True
                    print(f"Found trifluoromethyl sulfide group in molecule: {node['smiles']}")

                # Check for trifluorovinyl group
                trifluorovinyl_pattern = Chem.MolFromSmarts("[#6]=[C]([F])[F]")
                if mol.HasSubstructMatch(trifluorovinyl_pattern):
                    has_trifluorovinyl = True
                    print(f"Found trifluorovinyl group in molecule: {node['smiles']}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Traverse the route
    dfs_traverse(route)

    # Check if both fluorinated building blocks are present
    strategy_detected = has_trifluoromethyl_sulfide and has_trifluorovinyl

    if strategy_detected:
        print(
            "Detected fluorinated building blocks strategy with both trifluoromethyl sulfide and trifluorovinyl groups"
        )
    else:
        print("Did not detect complete fluorinated building blocks strategy")
        print(
            f"Trifluoromethyl sulfide: {has_trifluoromethyl_sulfide}, Trifluorovinyl: {has_trifluorovinyl}"
        )

    return strategy_detected
