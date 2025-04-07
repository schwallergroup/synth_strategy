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
    Detects if the synthetic route involves multiple nitrogen-containing functional groups
    (amines, amides, carbamates).
    """
    has_multiple_n_groups = False

    def dfs_traverse(node):
        nonlocal has_multiple_n_groups

        if (
            node["type"] == "mol"
            and "smiles" in node
            and not node.get("in_stock", False)
        ):
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if not mol:
                    return

                # Define patterns for different N-containing functional groups
                amine_pattern = Chem.MolFromSmarts("[NX3;H0,H1,H2]")
                amide_pattern = Chem.MolFromSmarts("[NX3][CX3]=[OX1]")
                carbamate_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[OX2]")

                # Count different N-containing functional groups
                n_groups = 0
                if mol.HasSubstructMatch(amine_pattern):
                    n_groups += 1
                if mol.HasSubstructMatch(amide_pattern):
                    n_groups += 1
                if mol.HasSubstructMatch(carbamate_pattern):
                    n_groups += 1

                if n_groups >= 2:
                    print(
                        f"Found molecule with {n_groups} different N-containing functional groups"
                    )
                    has_multiple_n_groups = True
            except:
                pass  # Handle parsing errors gracefully

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_multiple_n_groups
