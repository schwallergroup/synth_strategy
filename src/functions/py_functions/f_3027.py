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
    This function detects if the synthetic route uses O-benzylation for phenol protection.
    """
    benzylation_count = 0

    def dfs_traverse(node):
        nonlocal benzylation_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phenol in reactants
            phenol_pattern = Chem.MolFromSmarts("[OH]-[c]")
            benzyl_ether_pattern = Chem.MolFromSmarts("[O]-[CH2]-[c]")

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(phenol_pattern):
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(
                        benzyl_ether_pattern
                    ):
                        benzylation_count += 1
                        print(f"O-Benzylation detected: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Total O-benzylation reactions: {benzylation_count}")
    return benzylation_count >= 1
