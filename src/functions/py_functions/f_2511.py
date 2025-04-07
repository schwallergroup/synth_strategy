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
    This function detects a strategy involving silyl protection and deprotection.
    It looks for a silyl group that is present in intermediates but removed in the final product.
    """
    has_silyl_intermediate = False
    final_product_has_silyl = False

    def dfs_traverse(node, depth=0):
        nonlocal has_silyl_intermediate, final_product_has_silyl

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"]) if node["smiles"] else None

            if mol:
                # Pattern for silyl group
                silyl_pattern = Chem.MolFromSmarts("[Si]([C])([C])[C]")
                has_silyl = mol.HasSubstructMatch(silyl_pattern)

                if depth > 0:  # Intermediate
                    if has_silyl:
                        has_silyl_intermediate = True
                        print(
                            f"Silyl group detected in intermediate at depth {depth}: {node['smiles']}"
                        )
                else:  # Final product (depth 0)
                    if has_silyl:
                        final_product_has_silyl = True
                        print(
                            f"Silyl group detected in final product: {node['smiles']}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if silyl groups are present in intermediates but not in the final product
    return has_silyl_intermediate and not final_product_has_silyl
