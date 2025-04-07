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
    This function detects a strategy involving multiple carbonyl reduction steps
    in the synthesis route.
    """
    carbonyl_reduction_count = 0

    def dfs_traverse(node):
        nonlocal carbonyl_reduction_count

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a carbonyl reduction reaction
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            # Look for carbonyl pattern in reactant
                            carbonyl_pattern = Chem.MolFromSmarts("[C;!$(C=C)]=O")
                            if reactant_mol.HasSubstructMatch(carbonyl_pattern):
                                # Check if product has methyl or methylene where carbonyl was
                                methyl_pattern = Chem.MolFromSmarts("[C;!$(C=C)][C;H2,H3]")
                                if product_mol.HasSubstructMatch(methyl_pattern):
                                    carbonyl_reduction_count += 1
                                    print(
                                        f"Found carbonyl reduction at depth: {node.get('depth', 'unknown')}"
                                    )
                except Exception as e:
                    print(f"Error in carbonyl reduction detection: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Total carbonyl reductions found: {carbonyl_reduction_count}")
    return carbonyl_reduction_count >= 2
