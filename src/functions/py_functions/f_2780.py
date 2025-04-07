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
    This function detects a strategy involving multiple C-N bond formations.
    It counts the number of reactions that form or break C-N bonds.
    """
    c_n_bond_formations = 0

    def dfs_traverse(node):
        nonlocal c_n_bond_formations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                    # Define C-N bond pattern
                    c_n_bond_pattern = Chem.MolFromSmarts("[#6]-[#7]")

                    # Count C-N bonds in product and reactants
                    if product_mol:
                        product_c_n_bonds = len(
                            product_mol.GetSubstructMatches(c_n_bond_pattern)
                        )
                    else:
                        product_c_n_bonds = 0

                    reactant_c_n_bonds = sum(
                        len(r.GetSubstructMatches(c_n_bond_pattern))
                        for r in reactant_mols
                        if r
                    )

                    # In retrosynthesis, if product has more C-N bonds than reactants combined,
                    # it means a C-N bond was broken (formed in forward synthesis)
                    if product_c_n_bonds > reactant_c_n_bonds:
                        c_n_bond_formations += 1
                        print(
                            f"C-N bond formation detected, total count: {c_n_bond_formations}"
                        )
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if at least 2 C-N bond formations are found
    return c_n_bond_formations >= 2
