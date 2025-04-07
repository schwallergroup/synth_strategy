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
    This function detects a reductive amination strategy where a C-N bond
    is formed between a cyclic ketone and an amine.
    """
    reductive_amination_found = False

    def dfs_traverse(node, depth=0):
        nonlocal reductive_amination_found

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carbonyl pattern in reactants
            carbonyl_pattern = "C(=O)C"
            amine_pattern = "N([C,H])[C,H]"

            # Check for C-N bond in product where carbonyl was
            cn_bond_pattern = "CN([C,H])[C,H]"

            # Check if there's a carbonyl in reactants and a C-N bond in product
            has_carbonyl = any(
                Chem.MolFromSmiles(r)
                and Chem.MolFromSmiles(r).HasSubstructMatch(
                    Chem.MolFromSmarts(carbonyl_pattern)
                )
                for r in reactants
                if Chem.MolFromSmiles(r)
            )

            has_amine = any(
                Chem.MolFromSmiles(r)
                and Chem.MolFromSmiles(r).HasSubstructMatch(
                    Chem.MolFromSmarts(amine_pattern)
                )
                for r in reactants
                if Chem.MolFromSmiles(r)
            )

            product_mol = Chem.MolFromSmiles(product)
            has_cn_bond = product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts(cn_bond_pattern)
            )

            if has_carbonyl and has_amine and has_cn_bond:
                print(f"Found reductive amination at depth {depth}")
                reductive_amination_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return reductive_amination_found
