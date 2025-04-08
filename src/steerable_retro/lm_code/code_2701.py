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
    Detects if the route involves formation of an isocyanate from an amine.
    """
    isocyanate_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal isocyanate_formation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Isocyanate pattern
            isocyanate_pattern = Chem.MolFromSmarts("[#8]=[#6]=[#7]")

            # Amine pattern
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            if product is not None and isocyanate_pattern is not None:
                if product.HasSubstructMatch(isocyanate_pattern):
                    # Check if a reactant has an amine
                    for r in reactants:
                        if (
                            r is not None
                            and amine_pattern is not None
                            and r.HasSubstructMatch(amine_pattern)
                        ):
                            print("Detected isocyanate formation from amine at depth", depth)
                            isocyanate_formation_detected = True
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return isocyanate_formation_detected
