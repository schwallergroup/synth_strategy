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
    This function detects a strategy where halogenated aromatic systems (with F and Cl)
    are maintained throughout the synthesis.
    """
    has_halogenated_aromatic_target = False
    has_halogenated_aromatic_starting_material = False

    def dfs_traverse(node, depth=0):
        nonlocal has_halogenated_aromatic_target, has_halogenated_aromatic_starting_material

        if node["type"] == "mol":
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for halogenated aromatic pattern
                    f_aromatic_pattern = Chem.MolFromSmarts("c-[F]")
                    cl_aromatic_pattern = Chem.MolFromSmarts("c-[Cl]")

                    has_f = mol.HasSubstructMatch(f_aromatic_pattern)
                    has_cl = mol.HasSubstructMatch(cl_aromatic_pattern)

                    if has_f and has_cl:
                        if depth == 0:  # Target molecule
                            has_halogenated_aromatic_target = True
                            print(f"Target molecule has F and Cl substituted aromatic system")
                        elif node.get("in_stock", False):  # Starting material
                            has_halogenated_aromatic_starting_material = True
                            print(f"Starting material has F and Cl substituted aromatic system")
            except Exception as e:
                print(f"Error processing molecule: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both target and starting material have halogenated aromatics
    return has_halogenated_aromatic_target and has_halogenated_aromatic_starting_material
