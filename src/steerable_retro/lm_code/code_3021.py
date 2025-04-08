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
    Detects if the synthesis route involves nucleophilic substitution reactions,
    particularly focusing on C-N and C-S bond formations where a halogen is replaced.
    """
    nucleophilic_substitutions = 0

    def dfs_traverse(node, depth=0):
        nonlocal nucleophilic_substitutions

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and all(r for r in reactant_mols):
                # Check for C-N bond formation (where C was previously bonded to halogen)
                c_n_pattern = Chem.MolFromSmarts("[#6]-[#7]")
                c_x_pattern = Chem.MolFromSmarts("[#6]-[#9,#17,#35,#53]")  # C-halogen

                # Check for C-S bond formation (where C was previously bonded to halogen)
                c_s_pattern = Chem.MolFromSmarts("[#6]-[#16]")

                # Find C-N bonds in product that weren't in reactants
                if product_mol.HasSubstructMatch(c_n_pattern) and any(
                    r.HasSubstructMatch(c_x_pattern) for r in reactant_mols
                ):
                    nucleophilic_substitutions += 1
                    print(f"C-N nucleophilic substitution detected at depth {depth}")

                # Find C-S bonds in product that weren't in reactants
                if product_mol.HasSubstructMatch(c_s_pattern) and any(
                    r.HasSubstructMatch(c_x_pattern) for r in reactant_mols
                ):
                    nucleophilic_substitutions += 1
                    print(f"C-S nucleophilic substitution detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    has_multiple_substitutions = nucleophilic_substitutions >= 2
    if has_multiple_substitutions:
        print(f"Multiple nucleophilic substitutions detected: {nucleophilic_substitutions}")

    return has_multiple_substitutions
