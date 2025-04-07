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
    Detects a synthetic strategy involving reduction of a diester to a diol.
    """
    has_diester_reduction = False

    def dfs_traverse(node):
        nonlocal has_diester_reduction

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                if reactants and product:
                    # Check for diester pattern in reactants
                    diester_pattern = Chem.MolFromSmarts(
                        "[C](=[O])[O][CH2][CH3].[C](=[O])[O][CH2][CH3]"
                    )

                    # Check for diol pattern in product
                    diol_pattern = Chem.MolFromSmarts("[OH][CH2][CH][CH2][OH]")

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            diester_pattern
                        ):
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol and product_mol.HasSubstructMatch(
                                diol_pattern
                            ):
                                has_diester_reduction = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Diester reduction to diol strategy: {has_diester_reduction}")
    return has_diester_reduction
