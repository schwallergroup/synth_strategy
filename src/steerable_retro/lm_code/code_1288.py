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
    Detects a strategy where a nitro group is reduced to an amine and then
    further functionalized to a tertiary amine.
    """
    # Track if we found the required transformations
    found_nitro_to_amine = False
    found_amine_to_tertiary = False

    def dfs_traverse(node):
        nonlocal found_nitro_to_amine, found_amine_to_tertiary

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Parse molecules
                try:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                    product = Chem.MolFromSmiles(product_smiles)

                    # Check for nitro to amine transformation
                    nitro_pattern = Chem.MolFromSmarts("[#7+](=[O])[O-]")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")

                    reactants_have_nitro = any(
                        r.HasSubstructMatch(nitro_pattern) for r in reactants if r is not None
                    )
                    product_has_amine = product is not None and product.HasSubstructMatch(
                        amine_pattern
                    )

                    if reactants_have_nitro and product_has_amine:
                        print("Found nitro to amine transformation")
                        found_nitro_to_amine = True

                    # Check for amine to tertiary amine transformation
                    tertiary_amine_pattern = Chem.MolFromSmarts("[#7]([#6])([#6])")

                    reactants_have_amine = any(
                        r.HasSubstructMatch(amine_pattern) for r in reactants if r is not None
                    )
                    product_has_tertiary_amine = product is not None and product.HasSubstructMatch(
                        tertiary_amine_pattern
                    )

                    if reactants_have_amine and product_has_tertiary_amine:
                        print("Found amine to tertiary amine transformation")
                        found_amine_to_tertiary = True

                except:
                    print(f"Error processing reaction SMILES: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if both transformations are found
    return found_nitro_to_amine and found_amine_to_tertiary
