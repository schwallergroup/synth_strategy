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
    Detects if the synthesis includes sequential alkylation of a nitrogen atom,
    particularly in the context of a piperidine.
    """
    n_alkylation_count = 0
    sequential_alkylation_detected = False

    def dfs_traverse(node):
        nonlocal n_alkylation_count, sequential_alkylation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants = reactants_part.split(".")

                # Check for piperidine in reactants
                piperidine_pattern = Chem.MolFromSmarts("N1CCCCC1")

                # Check for alkylation reaction
                for reactant in reactants:
                    r_mol = Chem.MolFromSmiles(reactant)
                    if r_mol and r_mol.HasSubstructMatch(piperidine_pattern):
                        product_mol = Chem.MolFromSmiles(product_part)

                        # Check if the product has an alkylated nitrogen
                        alkylated_n_pattern = Chem.MolFromSmarts("N1(C*)CCCCC1")
                        if product_mol and product_mol.HasSubstructMatch(
                            alkylated_n_pattern
                        ):
                            n_alkylation_count += 1
                            print(
                                f"Detected nitrogen alkylation step (count: {n_alkylation_count})"
                            )

                            if n_alkylation_count >= 2:
                                sequential_alkylation_detected = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return sequential_alkylation_detected
