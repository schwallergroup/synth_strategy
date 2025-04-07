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
    Detects if the synthesis route involves reduction of a nitro group to an amine.
    """
    nitro_reduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro group in reactants
            nitro_pattern = Chem.MolFromSmarts("[#6]-[NX3+](=[OX1])[OX1-]")

            # Check for amine in product where nitro was
            amine_pattern = Chem.MolFromSmarts("[NX3]")

            product_mol = Chem.MolFromSmiles(product)

            # Check if reactants contain nitro and product contains new amine
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(nitro_pattern):
                    if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                        # This is a simplification - ideally we'd check that the amine is at the same position
                        # where the nitro group was, but that requires more complex analysis
                        nitro_reduction_found = True
                        print(f"Nitro reduction detected at depth {depth}")
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return nitro_reduction_found
