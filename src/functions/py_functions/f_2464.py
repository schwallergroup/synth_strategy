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
    Detects ester deprotection to form amide or carboxylic acid, a common protecting group strategy.
    """
    ester_deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal ester_deprotection_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for methyl ester pattern in reactants
                methyl_ester_pattern = Chem.MolFromSmarts("[#6](=[#8])[O][C]")

                # Check for carboxylic acid or amide pattern in product
                acid_pattern = Chem.MolFromSmarts("[#6](=[#8])[O][H]")
                amide_pattern = Chem.MolFromSmarts("[#6](=[#8])[N]")

                has_methyl_ester = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(methyl_ester_pattern)
                    for r in reactants
                )

                product_mol = Chem.MolFromSmiles(product)
                has_acid_or_amide = product_mol is not None and (
                    product_mol.HasSubstructMatch(acid_pattern)
                    or product_mol.HasSubstructMatch(amide_pattern)
                )

                if has_methyl_ester and has_acid_or_amide:
                    ester_deprotection_found = True
                    print(f"Found ester deprotection at depth {depth}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return ester_deprotection_found
