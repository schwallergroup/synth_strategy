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
    This function detects a strategy involving multiple C-N bond formations throughout the synthesis.
    """
    cn_bond_formations = 0

    def dfs_traverse(node):
        nonlocal cn_bond_formations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
                product_mol = Chem.MolFromSmiles(product_part)

                if product_mol:
                    # Check for C-N bonds in product that aren't in reactants
                    cn_bond_pattern = Chem.MolFromSmarts("[#6]-[#7]")
                    product_cn_matches = product_mol.GetSubstructMatches(cn_bond_pattern)

                    # Count C-N bonds in reactants
                    reactant_cn_count = 0
                    for r_mol in reactants_mols:
                        if r_mol:
                            reactant_cn_count += len(r_mol.GetSubstructMatches(cn_bond_pattern))

                    # If product has more C-N bonds than reactants combined, a C-N bond was formed
                    if len(product_cn_matches) > reactant_cn_count:
                        cn_bond_formations += 1
                        print(f"Found C-N bond formation")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If at least 2 C-N bond formations are detected
    return cn_bond_formations >= 2
