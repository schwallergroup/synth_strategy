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
    This function detects if the synthesis route involves alkylation of a phenol.
    """
    phenol_alkylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal phenol_alkylation_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phenol alkylation: phenol to alkoxybenzene
            phenol_pattern = Chem.MolFromSmarts("[OH]-[c]")
            alkoxybenzene_pattern = Chem.MolFromSmarts("[c]-[O]-[CH2]")

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(alkoxybenzene_pattern):
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(phenol_pattern):
                        print(f"Phenol alkylation detected at depth {depth}")
                        phenol_alkylation_found = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return phenol_alkylation_found
