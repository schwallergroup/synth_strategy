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
    This function detects if the synthetic route contains an N-alkylation reaction
    of a heterocyclic NH with an alkyl halide.
    """
    n_alkylation_detected = False

    def dfs_traverse(node):
        nonlocal n_alkylation_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for heterocyclic NH and alkyl halide in reactants
            nh_pattern = Chem.MolFromSmarts("[nH]")
            alkyl_halide_pattern = Chem.MolFromSmarts("[#6]-[#35,#53,#17]")

            # Check for N-alkylated product
            n_alkylated_pattern = Chem.MolFromSmarts("[#6]-[n]")

            has_nh = False
            has_alkyl_halide = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(nh_pattern):
                        has_nh = True
                    if mol.HasSubstructMatch(alkyl_halide_pattern):
                        has_alkyl_halide = True

            product_mol = Chem.MolFromSmiles(product)
            has_n_alkylated = product_mol and product_mol.HasSubstructMatch(
                n_alkylated_pattern
            )

            if has_nh and has_alkyl_halide and has_n_alkylated:
                print("N-alkylation detected")
                n_alkylation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return n_alkylation_detected
