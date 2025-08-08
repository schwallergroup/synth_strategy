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
    This function detects the construction of a complex heterocyclic system
    with multiple heterocycles (e.g., thiazole and pyridine).
    """
    # Track heterocycles at each depth
    heterocycles_by_depth = {}

    def dfs_traverse(node):
        if node["type"] == "mol":
            smiles = node["smiles"]
            depth = node.get("depth", -1)

            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Define patterns for different heterocycles
                thiazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#16][#6][#6]1")
                pyridine_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#7][#6]1")

                # Count heterocycles
                thiazole_count = len(mol.GetSubstructMatches(thiazole_pattern))
                pyridine_count = len(mol.GetSubstructMatches(pyridine_pattern))

                heterocycles_by_depth[depth] = {
                    "thiazole": thiazole_count,
                    "pyridine": pyridine_count,
                    "total": thiazole_count + pyridine_count,
                }

                if thiazole_count > 0 or pyridine_count > 0:
                    print(
                        f"At depth {depth}: Found {thiazole_count} thiazoles and {pyridine_count} pyridines"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Analyze if we're constructing a complex heterocyclic system
    depths = sorted(heterocycles_by_depth.keys())

    if not depths:
        return False

    # Check if the final product has multiple heterocycles
    min_depth = min(depths)
    complex_system = heterocycles_by_depth[min_depth]["total"] >= 2

    # Check if heterocycles are being added during synthesis
    if len(depths) > 1:
        max_depth = max(depths)
        heterocycle_construction = (
            heterocycles_by_depth[min_depth]["total"] > heterocycles_by_depth[max_depth]["total"]
        )

        if heterocycle_construction:
            print("Heterocycles are being constructed during synthesis")

    print(f"Complex heterocyclic system construction strategy detected: {complex_system}")
    return complex_system
