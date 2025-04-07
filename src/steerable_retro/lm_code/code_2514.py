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
    This function detects if the synthesis involves O-alkylation of a phenol.
    """
    o_alkylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal o_alkylation_detected

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    prod_mol = Chem.MolFromSmiles(product)

                    # Check for O-alkylation of phenol
                    phenol_patt = Chem.MolFromSmarts("[OH]c1ccccc1")
                    alkylated_phenol_patt = Chem.MolFromSmarts("[O]([CH2][CH])c1ccccc1")

                    for reactant in reactants:
                        react_mol = Chem.MolFromSmiles(reactant)
                        if react_mol and prod_mol:
                            if react_mol.HasSubstructMatch(
                                phenol_patt
                            ) and prod_mol.HasSubstructMatch(alkylated_phenol_patt):
                                print(f"O-alkylation of phenol detected at depth {depth}")
                                o_alkylation_detected = True
                except Exception as e:
                    print(f"Error in SMILES processing: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return o_alkylation_detected
