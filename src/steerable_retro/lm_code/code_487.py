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
    This function detects if the synthesis involves amide bond formation using an acid chloride.
    """
    amide_formation_found = False

    # SMARTS patterns
    acid_chloride_pattern = Chem.MolFromSmarts("[C](=[O])-[Cl]")
    amine_pattern = Chem.MolFromSmarts("[N;!$(N=*);!$(NC=O)]")
    amide_pattern = Chem.MolFromSmarts("[C](=[O])-[N]")

    def dfs_traverse(node):
        nonlocal amide_formation_found

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
            product = Chem.MolFromSmiles(products_part)

            # Check for amide formation from acid chloride and amine
            if (
                any(r is not None and r.HasSubstructMatch(acid_chloride_pattern) for r in reactants)
                and any(r is not None and r.HasSubstructMatch(amine_pattern) for r in reactants)
                and product is not None
                and product.HasSubstructMatch(amide_pattern)
            ):
                amide_formation_found = True
                print("Found amide formation from acid chloride")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return amide_formation_found
