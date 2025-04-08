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
    Detects if the route involves multiple C-N bond formation steps
    """
    cn_formations = 0

    def dfs_traverse(node):
        nonlocal cn_formations

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactants_mol = Chem.MolFromSmiles(reactants)
                product_mol = Chem.MolFromSmiles(product)

                if reactants_mol and product_mol:
                    # Check for C-N bond formation by comparing product to reactants
                    cn_pattern = Chem.MolFromSmarts("[#6]-[#7]")

                    reactant_cn_count = (
                        len(reactants_mol.GetSubstructMatches(cn_pattern))
                        if reactants_mol.HasSubstructMatch(cn_pattern)
                        else 0
                    )
                    product_cn_count = (
                        len(product_mol.GetSubstructMatches(cn_pattern))
                        if product_mol.HasSubstructMatch(cn_pattern)
                        else 0
                    )

                    if product_cn_count > reactant_cn_count:
                        print(f"Found C-N bond formation: {reactant_cn_count} → {product_cn_count}")
                        cn_formations += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return cn_formations >= 2  # Return True if at least 2 C-N bond formations
