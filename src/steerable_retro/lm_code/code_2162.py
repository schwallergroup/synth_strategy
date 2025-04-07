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
    Detects if the synthesis route involves multiple Boc deprotection steps.
    """
    boc_deprotection_count = 0

    def dfs_traverse(node):
        nonlocal boc_deprotection_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactant has Boc group and product doesn't
            boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)[N]")
            nh_pattern = Chem.MolFromSmarts("[NH]")

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if not reactant_mol:
                    continue

                product_mol = Chem.MolFromSmiles(product)
                if not product_mol:
                    continue

                if (
                    reactant_mol.HasSubstructMatch(boc_pattern)
                    and not product_mol.HasSubstructMatch(boc_pattern)
                    and product_mol.HasSubstructMatch(nh_pattern)
                ):
                    boc_deprotection_count += 1
                    print(f"Found Boc deprotection at depth: {node.get('depth', 'unknown')}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return boc_deprotection_count >= 2
