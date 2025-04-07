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
    This function detects a synthetic strategy involving multiple protection/deprotection steps.
    """
    protection_count = 0

    def dfs_traverse(node):
        nonlocal protection_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc deprotection
            boc_pattern = Chem.MolFromSmarts("[#6]OC(=O)N")
            product_mol = Chem.MolFromSmiles(product)

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if (
                    reactant_mol
                    and boc_pattern
                    and reactant_mol.HasSubstructMatch(boc_pattern)
                ):
                    if product_mol and not product_mol.HasSubstructMatch(boc_pattern):
                        print("Found Boc deprotection")
                        protection_count += 1

            # Check for benzoate deprotection
            benzoate_pattern = Chem.MolFromSmarts("c1ccccc1C(=O)O[#6]")
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if (
                    reactant_mol
                    and benzoate_pattern
                    and reactant_mol.HasSubstructMatch(benzoate_pattern)
                ):
                    if product_mol and not product_mol.HasSubstructMatch(
                        benzoate_pattern
                    ):
                        print("Found benzoate deprotection")
                        protection_count += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return protection_count >= 2  # At least two protection/deprotection steps
