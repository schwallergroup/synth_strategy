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
    This function detects if the synthetic route involves building nitrogen-rich
    heterocyclic scaffolds (containing 3+ nitrogen atoms).
    """
    nitrogen_rich_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrogen_rich_detected

        if node["type"] == "mol" and not node.get("in_stock", False):
            # Count nitrogen atoms in the molecule
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Count nitrogen atoms
                nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

                # Check for heterocycles
                ring_info = mol.GetRingInfo()
                has_nitrogen_ring = False

                for ring in ring_info.AtomRings():
                    ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
                    if any(atom.GetAtomicNum() == 7 for atom in ring_atoms):
                        has_nitrogen_ring = True
                        break

                if nitrogen_count >= 3 and has_nitrogen_ring:
                    print(f"Detected nitrogen-rich heterocycle at depth {depth}")
                    nitrogen_rich_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if nitrogen_rich_detected:
        print("Nitrogen-rich heterocycle strategy detected")
    return nitrogen_rich_detected
