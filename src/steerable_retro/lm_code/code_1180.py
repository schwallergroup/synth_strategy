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
    This function detects if the synthetic route involves silyl protection of a hydroxyl group.
    """
    silyl_protection_found = False

    def dfs_traverse(node):
        nonlocal silyl_protection_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactants contain a hydroxyl group and product contains a silyl ether
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol:
                # Check for hydroxyl in reactants
                hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
                has_hydroxyl = any(
                    mol.HasSubstructMatch(hydroxyl_pattern) for mol in reactant_mols if mol
                )

                # Check for silyl ether in product
                silyl_pattern = Chem.MolFromSmarts("[OX2][Si]")
                has_silyl = product_mol.HasSubstructMatch(silyl_pattern) if product_mol else False

                if has_hydroxyl and has_silyl:
                    print("Silyl protection detected")
                    silyl_protection_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return silyl_protection_found
