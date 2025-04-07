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
    This function detects if the synthesis includes a nitro reduction to amine.
    """
    nitro_reduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a nitro reduction reaction
                nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                # Check if nitro group in reactant and amine in product
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    product_mol = Chem.MolFromSmiles(product)

                    if (
                        reactant_mol
                        and product_mol
                        and reactant_mol.HasSubstructMatch(nitro_pattern)
                        and product_mol.HasSubstructMatch(amine_pattern)
                        and not product_mol.HasSubstructMatch(nitro_pattern)
                    ):
                        nitro_reduction_found = True
                        print(f"Nitro reduction found at depth {depth}")
                        break

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return nitro_reduction_found
