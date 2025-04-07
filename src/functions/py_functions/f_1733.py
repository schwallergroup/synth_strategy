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
    Detects if the synthesis route involves N-methylation of a heterocyclic system.
    """
    n_methylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal n_methylation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for N-methylation pattern
            product_mol = Chem.MolFromSmiles(product)
            n_methyl_pattern = Chem.MolFromSmarts("[N;R][C;H3]")

            if product_mol and product_mol.HasSubstructMatch(n_methyl_pattern):
                # Check if this is a methylation (methyl addition)
                has_n_methyl_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        n_methyl_pattern
                    ):
                        has_n_methyl_in_reactants = True
                        break

                if not has_n_methyl_in_reactants:
                    n_methylation_detected = True
                    print(f"N-methylation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return n_methylation_detected
