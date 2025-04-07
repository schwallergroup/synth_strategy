#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects if the final product contains multiple (>3) connected aromatic rings.
    """
    has_multi_aromatic = False

    def dfs_traverse(node, is_root=True):
        nonlocal has_multi_aromatic

        if node["type"] == "mol" and is_root:  # Final product node (root)
            try:
                mol_smiles = node["smiles"]
                print(f"Analyzing final product: {mol_smiles}")
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Get all aromatic rings
                    aromatic_rings = []
                    ring_info = mol.GetRingInfo()
                    for ring_atoms in ring_info.AtomRings():
                        is_aromatic = True
                        for atom_idx in ring_atoms:
                            atom = mol.GetAtomWithIdx(atom_idx)
                            if not atom.GetIsAromatic():
                                is_aromatic = False
                                break
                        if is_aromatic:
                            aromatic_rings.append(set(ring_atoms))

                    print(f"Found {len(aromatic_rings)} aromatic rings")

                    # Check if we have more than 3 aromatic rings
                    if len(aromatic_rings) > 3:
                        print(
                            f"Final product contains more than 3 aromatic rings: {mol_smiles}"
                        )
                        has_multi_aromatic = True
                        return

                    # Check if rings are connected (share atoms or bonds)
                    if len(aromatic_rings) >= 2:
                        # Build a graph of connected rings
                        connected_rings = {}
                        for i in range(len(aromatic_rings)):
                            connected_rings[i] = []
                            for j in range(len(aromatic_rings)):
                                if i != j and aromatic_rings[i].intersection(
                                    aromatic_rings[j]
                                ):
                                    connected_rings[i].append(j)

                        # Count connected components and their sizes
                        visited = set()
                        connected_components = []

                        def dfs_rings(ring_idx, component):
                            visited.add(ring_idx)
                            component.append(ring_idx)
                            for neighbor in connected_rings[ring_idx]:
                                if neighbor not in visited:
                                    dfs_rings(neighbor, component)

                        for i in range(len(aromatic_rings)):
                            if i not in visited:
                                component = []
                                dfs_rings(i, component)
                                connected_components.append(component)

                        # Check if any component has more than 3 rings
                        for component in connected_components:
                            print(
                                f"Found connected component with {len(component)} rings"
                            )
                            if len(component) > 3:
                                print(
                                    f"Final product contains a connected system of more than 3 aromatic rings"
                                )
                                has_multi_aromatic = True
                                return
            except Exception as e:
                print(f"Error analyzing molecule: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, is_root=False)

    dfs_traverse(route)
    print(f"Multi-aromatic system found: {has_multi_aromatic}")
    return has_multi_aromatic
