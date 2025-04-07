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
    Detects if the synthesis route involves the formation of a heterocyclic ring,
    particularly a triazole ring.
    """
    # Track if we found a heterocycle formation reaction
    found_heterocycle_formation = False
    formation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle_formation, formation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this reaction forms a triazole ring
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol:
                # Check for triazole pattern in product
                triazole_pattern = Chem.MolFromSmarts("[n]1[n][c][c][n]1")
                if product_mol.HasSubstructMatch(triazole_pattern):
                    # Check if reactants don't have triazole
                    has_triazole_in_reactants = False
                    for mol in reactant_mols:
                        if mol and mol.HasSubstructMatch(triazole_pattern):
                            has_triazole_in_reactants = True
                            break

                    if not has_triazole_in_reactants:
                        found_heterocycle_formation = True
                        formation_depth = depth
                        print(f"Found triazole formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    if found_heterocycle_formation:
        print(
            f"Heterocycle formation strategy detected: Triazole formation at depth {formation_depth}"
        )
    else:
        print("Heterocycle formation strategy not detected")

    return found_heterocycle_formation
