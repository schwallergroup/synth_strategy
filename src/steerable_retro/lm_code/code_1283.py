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
    Detects a synthetic strategy involving biaryl formation via cross-coupling.
    """
    has_biaryl_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_biaryl_formation

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for biaryl formation (two aromatic rings connected)
            # Look for reactions where separate aromatic rings in reactants become connected in product
            if len(reactants) >= 2:
                # Check if product contains biaryl motif
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol:
                        biaryl_pattern = Chem.MolFromSmarts(
                            "[c]!@[c]"
                        )  # Non-fused aromatic C connected to another aromatic C
                        if prod_mol.HasSubstructMatch(biaryl_pattern):
                            # Check if reactants have aromatic rings
                            aromatic_reactants = 0
                            for r in reactants:
                                r_mol = Chem.MolFromSmiles(r)
                                if r_mol and r_mol.HasSubstructMatch(Chem.MolFromSmarts("c")):
                                    aromatic_reactants += 1

                            if aromatic_reactants >= 2:
                                has_biaryl_formation = True
                                print(f"Found biaryl formation at depth {depth}")
                except:
                    pass  # Handle parsing errors gracefully

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_biaryl_formation
