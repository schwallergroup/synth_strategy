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
    Detects if the synthetic route uses iodine as a persistent functional handle
    while fluorine is used as a leaving group.
    """
    iodine_persistent = False
    fluorine_leaving_group = False

    def dfs_traverse(node):
        nonlocal iodine_persistent, fluorine_leaving_group

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactants_mol = Chem.MolFromSmiles(".".join(reactants))
            product_mol = Chem.MolFromSmiles(product)

            if reactants_mol and product_mol:
                # Check for iodine persistence
                iodine_pattern = Chem.MolFromSmarts("[#53]")
                if reactants_mol.HasSubstructMatch(
                    iodine_pattern
                ) and product_mol.HasSubstructMatch(iodine_pattern):
                    iodine_persistent = True

                # Check for fluorine as leaving group
                fluoro_aromatic_pattern = Chem.MolFromSmarts("c[F]")
                if len(
                    reactants_mol.GetSubstructMatches(fluoro_aromatic_pattern)
                ) > len(product_mol.GetSubstructMatches(fluoro_aromatic_pattern)):
                    fluorine_leaving_group = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if both conditions are met
    result = iodine_persistent and fluorine_leaving_group
    print(f"Halogen-directed synthesis strategy detected: {result}")
    return result
