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
    This function detects if the route involves an adamantane-like polycyclic scaffold.
    """
    adamantane_present = False

    def is_adamantane_like(smiles):
        """
        Helper function to detect adamantane-like polycyclic scaffolds.
        Adamantane is a tricyclic molecule with a cage-like structure.
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False

            # Check for tricyclic or tetracyclic structure with specific ring patterns
            # Adamantane has 3 fused cyclohexane rings in a specific arrangement
            ri = mol.GetRingInfo()

            # Check if molecule has at least 3 rings
            if ri.NumRings() < 3:
                return False

            # Check for cyclohexane rings
            cyclohexane_count = 0
            for ring_atoms in ri.AtomRings():
                if len(ring_atoms) == 6:  # Cyclohexane rings have 6 atoms
                    cyclohexane_count += 1

            # Adamantane-like structures typically have at least 3 cyclohexane rings
            if cyclohexane_count >= 3:
                # Check for cage-like structure by examining the connectivity
                # In adamantane, each carbon is connected to 3 or 4 other carbons
                carbon_connectivity = []
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() == "C":
                        carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == "C"]
                        carbon_connectivity.append(len(carbon_neighbors))

                # Count carbons with 3 or 4 carbon neighbors (characteristic of adamantane)
                highly_connected_carbons = sum(1 for conn in carbon_connectivity if conn >= 3)

                # Adamantane has several highly connected carbons
                if highly_connected_carbons >= 4:
                    print(f"Adamantane-like scaffold detected in molecule {smiles}")
                    return True

            return False
        except Exception as e:
            print(f"Error in is_adamantane_like: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal adamantane_present

        if node["type"] == "mol" and "smiles" in node:
            # Use the helper function to detect adamantane-like structure
            if is_adamantane_like(node["smiles"]):
                adamantane_present = True
                print(
                    f"Adamantane-like scaffold detected in molecule {node['smiles']} at depth {depth}"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return adamantane_present
