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
    This function detects a synthetic strategy involving multiple C-N bond formations
    throughout the synthesis (at least 2).
    """
    c_n_bond_formations = 0

    def dfs_traverse(node, depth=0):
        nonlocal c_n_bond_formations

        if node["type"] == "reaction":
            # Check if this reaction forms a C-N bond
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Convert to molecules
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                    ]

                    # Check for C-N bond formation
                    c_n_pattern = Chem.MolFromSmarts("[#6]-[#7]")
                    if product_mol:
                        product_c_n_bonds = len(product_mol.GetSubstructMatches(c_n_pattern))

                        # Count C-N bonds in reactants
                        reactant_c_n_bonds = sum(
                            len(mol.GetSubstructMatches(c_n_pattern))
                            for mol in reactant_mols
                            if mol
                        )

                        if product_c_n_bonds > reactant_c_n_bonds:
                            c_n_bond_formations += 1
                            print(f"C-N bond formation detected in reaction at depth {depth}")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if at least 2 C-N bond formations were detected
    return c_n_bond_formations >= 2
