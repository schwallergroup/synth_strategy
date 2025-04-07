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
    Detects if the synthesis involves an amino acid derivative with a benzyl side chain.
    """
    has_amino_acid_derivative = False

    def dfs_traverse(node):
        nonlocal has_amino_acid_derivative

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])

            if mol:
                # Pattern for amino acid with benzyl side chain
                amino_acid_pattern = Chem.MolFromSmarts(
                    "[NH][C]([C](=[O])[O,OC])[CH2]c1ccccc1"
                )

                if mol.HasSubstructMatch(amino_acid_pattern):
                    print("Found amino acid derivative with benzyl side chain")
                    has_amino_acid_derivative = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_amino_acid_derivative
