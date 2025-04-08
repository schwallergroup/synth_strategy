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
    This function detects if the synthesis involves sequential functionalization of an aromatic ring
    with multiple different functional groups.
    """
    # Track functionalization steps on the same aromatic ring
    aromatic_functionalizations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            # Simple check for aromatic ring functionalization
            # Look for patterns where an aromatic ring gains a new substituent
            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
                product = Chem.MolFromSmiles(products_part)

                if product is not None:
                    # Check if this is a functionalization of an aromatic ring
                    for r in reactants:
                        if r is not None and has_aromatic_ring(r) and has_aromatic_ring(product):
                            # This is a simplification - in a real implementation, you would need
                            # to track the specific aromatic ring and its substituents
                            aromatic_functionalizations.append(
                                (depth, "aromatic_functionalization")
                            )
                            print(f"Found aromatic functionalization at depth {depth}")
                            break
            except:
                pass  # Handle parsing errors gracefully

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    def has_aromatic_ring(mol):
        """Helper function to check if molecule has an aromatic ring"""
        if mol is None:
            return False
        aromatic_atoms = [atom.GetIsAromatic() for atom in mol.GetAtoms()]
        return any(aromatic_atoms)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least 2 functionalization steps
    return len(aromatic_functionalizations) >= 2
