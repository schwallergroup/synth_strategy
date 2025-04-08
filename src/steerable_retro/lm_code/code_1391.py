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
    Detects if the synthesis route involves reduction of an ester to an alcohol.
    """
    ester_reduction_found = False

    def dfs_traverse(node):
        nonlocal ester_reduction_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for ester reduction pattern
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactant_mol and product_mol:
                    # Check if reactant contains ester
                    ester_pattern = Chem.MolFromSmarts("[C](=[O])[O]")
                    # Check if product contains alcohol
                    alcohol_pattern = Chem.MolFromSmarts("[C]-[O;!H0]")

                    if (
                        reactant_mol.HasSubstructMatch(ester_pattern)
                        and product_mol.HasSubstructMatch(alcohol_pattern)
                        and not product_mol.HasSubstructMatch(ester_pattern)
                    ):
                        ester_reduction_found = True
                        print("Ester reduction to alcohol detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return ester_reduction_found
