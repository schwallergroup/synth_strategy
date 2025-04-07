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
    This function detects if the synthesis involves multiple C-O bond formations.
    """
    co_bond_formations = 0

    def dfs_traverse(node):
        nonlocal co_bond_formations

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check for C-O bond formation
            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
            product_mol = Chem.MolFromSmiles(product_part)

            if product_mol and all(r for r in reactants_mols):
                # Count C-O bonds in reactants and product
                co_pattern = Chem.MolFromSmarts("[#6]-[#8]")

                reactant_co_count = sum(
                    len(r.GetSubstructMatches(co_pattern)) for r in reactants_mols if r
                )
                product_co_count = len(product_mol.GetSubstructMatches(co_pattern))

                if product_co_count > reactant_co_count:
                    co_bond_formations += 1
                    print(
                        f"Found C-O bond formation: {reactant_co_count} -> {product_co_count}"
                    )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return co_bond_formations >= 2  # Return True if at least 2 C-O bond formations
