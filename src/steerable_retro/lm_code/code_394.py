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
    Detects a strategy involving protection of carboxylic acid as tert-butyl ester.
    """
    # Track if we found the protection
    protection_found = False

    # tert-butyl ester pattern
    tert_butyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC(C)(C)C")
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")

    def dfs_traverse(node, depth=0):
        nonlocal protection_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Convert to RDKit molecules
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                # Check for carboxylic acid protection
                if (
                    product is not None
                    and product.HasSubstructMatch(tert_butyl_ester_pattern)
                    and any(
                        r is not None and r.HasSubstructMatch(carboxylic_acid_pattern)
                        for r in reactants
                    )
                ):
                    protection_found = True
                    print(
                        f"Carboxylic acid protection as tert-butyl ester detected at depth {depth}"
                    )

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if protection_found:
        print("Carboxylic acid protection strategy detected")

    return protection_found
