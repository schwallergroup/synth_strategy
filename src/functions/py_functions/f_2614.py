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
    Detects if the synthesis route involves N-alkylation of a heterocycle
    (particularly pyrazole).
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

            # Check for N-alkylation pattern
            product_mol = Chem.MolFromSmiles(product)
            alkylated_n_patt = Chem.MolFromSmarts("[n;$(n1cccn1)][C;!$(C=O)]")

            if product_mol and product_mol.HasSubstructMatch(alkylated_n_patt):
                # Check if reactant had an unalkylated N
                unalkylated_n_patt = Chem.MolFromSmarts("[nH;$(n1cccn1)]")
                alkyl_halide_patt = Chem.MolFromSmarts("[C][Br,I,Cl]")

                has_unalkylated_n = False
                has_alkyl_halide = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    if reactant_mol.HasSubstructMatch(unalkylated_n_patt):
                        has_unalkylated_n = True
                    if reactant_mol.HasSubstructMatch(alkyl_halide_patt):
                        has_alkyl_halide = True

                if has_unalkylated_n and has_alkyl_halide:
                    print("N-alkylation of heterocycle detected")
                    n_alkylation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return n_alkylation_detected
