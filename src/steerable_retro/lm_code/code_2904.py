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
    Detects if the synthetic route involves late-stage N-alkylation of a heterocycle
    with a piperidine-containing fragment.
    """
    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction" and depth <= 1:  # Late stage (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains piperidine
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    piperidine_pattern = Chem.MolFromSmarts("[N]1[C][C][C][C][C]1")
                    if product_mol.HasSubstructMatch(piperidine_pattern):
                        # Check if this is an N-alkylation of a heterocycle
                        n_alkylation = False
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                # Look for aromatic nitrogen that's not alkylated
                                aromatic_n_pattern = Chem.MolFromSmarts("[n]")
                                if reactant_mol.HasSubstructMatch(aromatic_n_pattern):
                                    n_alkylation = True
                                    break

                        if n_alkylation:
                            print("Found late-stage N-alkylation with piperidine")
                            found_pattern = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_pattern
