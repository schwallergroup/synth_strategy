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
    This function detects late-stage nitro group reduction (depth 0-1).
    """
    late_nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal late_nitro_reduction

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if this is a nitro reduction
            if reactants_part and product_part:
                reactant_mol = Chem.MolFromSmiles(reactants_part)
                product_mol = Chem.MolFromSmiles(product_part)

                if reactant_mol and product_mol:
                    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")

                    if reactant_mol.HasSubstructMatch(
                        nitro_pattern
                    ) and product_mol.HasSubstructMatch(amine_pattern):
                        # Check if this is late stage (depth 0-1)
                        if depth <= 1:
                            late_nitro_reduction = True
                            print(
                                f"Detected late-stage nitro reduction at depth {depth}"
                            )

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return late_nitro_reduction
