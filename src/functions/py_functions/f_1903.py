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
    This function detects a late-stage fragment coupling strategy via C-O bond formation.
    It looks for a reaction in the first two steps that joins two complex fragments through an ether bond.
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Focus on late-stage reactions (depth 0 or 1)
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Check if we have multiple reactants (fragment coupling)
                reactants = reactants_part.split(".")
                if len(reactants) >= 2:
                    # Check for C-O bond formation (ether linkage)
                    product_mol = Chem.MolFromSmiles(product_part)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                    # Look for O atoms in the product that might be part of new ether bonds
                    for atom in product_mol.GetAtoms():
                        if atom.GetSymbol() == "O" and atom.GetDegree() == 2:
                            neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                            if (
                                len(neighbors) == 2
                            ):  # Ensure it's connected to two atoms
                                result = True
                                print(
                                    f"Found late-stage fragment coupling via C-O bond at depth {depth}"
                                )
                                return
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result
