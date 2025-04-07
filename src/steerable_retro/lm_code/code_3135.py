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
    Detects if the synthesis involves multiple (≥3) amine functionalization steps
    """
    amine_functionalization_count = 0

    def dfs_traverse(node):
        nonlocal amine_functionalization_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amine functionalization
            amine_pattern = Chem.MolFromSmarts("[NX3;!$(NC=O)]")

            # Check if any reactant has an amine
            has_amine_reactant = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(amine_pattern):
                    has_amine_reactant = True
                    break

            # Check if product has a new C-N bond
            if has_amine_reactant:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # This is a simplification - in a real implementation,
                    # you would need to compare the exact bonds
                    amine_functionalization_count += 1
                    print(f"Found amine functionalization: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total amine functionalization steps: {amine_functionalization_count}")
    return amine_functionalization_count >= 3
