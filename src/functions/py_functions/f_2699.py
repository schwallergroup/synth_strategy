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
    Detects if the route involves aromatic halogenation (specifically chlorination).
    """
    aromatic_halogenation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_halogenation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Aromatic C-Cl pattern
            aromatic_cl_pattern = Chem.MolFromSmarts("c[Cl]")

            if product is not None and aromatic_cl_pattern is not None:
                if product.HasSubstructMatch(aromatic_cl_pattern):
                    # Check if this bond is newly formed
                    reactants_have_pattern = False
                    for r in reactants:
                        if r is not None and r.HasSubstructMatch(aromatic_cl_pattern):
                            reactants_have_pattern = True
                            break

                    if not reactants_have_pattern:
                        print("Detected aromatic chlorination at depth", depth)
                        aromatic_halogenation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return aromatic_halogenation_detected
