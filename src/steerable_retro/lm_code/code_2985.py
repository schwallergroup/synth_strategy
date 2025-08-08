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
    Detects if the synthesis route involves compounds containing a spirocyclic system.
    A spirocyclic system consists of two rings that share exactly one atom (the spiro center).
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            try:
                # General detection of spirocyclic systems
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Get ring information
                    ring_info = mol.GetRingInfo()
                    if ring_info.NumRings() >= 2:  # Need at least 2 rings for a spirocyclic system
                        # Check all pairs of rings to see if they share exactly one atom
                        for i in range(ring_info.NumRings()):
                            ring_i_atoms = set(ring_info.AtomMembers(i))
                            for j in range(i + 1, ring_info.NumRings()):
                                ring_j_atoms = set(ring_info.AtomMembers(j))
                                # If rings share exactly one atom, it's a spiro system
                                shared_atoms = ring_i_atoms.intersection(ring_j_atoms)
                                if len(shared_atoms) == 1:
                                    print(
                                        f"Detected spirocyclic system at depth {depth}, SMILES: {mol_smiles}"
                                    )
                                    result = True
                                    return
            except Exception as e:
                print(f"Error processing molecule at depth {depth}: {e}")

        # Check reaction nodes to see if they create spirocyclic systems
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check if the product contains a spirocyclic system
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol:
                    ring_info = prod_mol.GetRingInfo()
                    if ring_info.NumRings() >= 2:
                        for i in range(ring_info.NumRings()):
                            ring_i_atoms = set(ring_info.AtomMembers(i))
                            for j in range(i + 1, ring_info.NumRings()):
                                ring_j_atoms = set(ring_info.AtomMembers(j))
                                shared_atoms = ring_i_atoms.intersection(ring_j_atoms)
                                if len(shared_atoms) == 1:
                                    print(
                                        f"Detected spirocyclic system in reaction product at depth {depth}, SMILES: {product}"
                                    )
                                    result = True
                                    return
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result
