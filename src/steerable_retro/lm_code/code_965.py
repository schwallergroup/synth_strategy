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
    This function detects multiple O-alkylation steps in the synthesis.
    """
    o_alkylation_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal o_alkylation_count

        if node["type"] == "reaction":
            # Get reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check for O-alkylation pattern
            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
            product_mol = Chem.MolFromSmiles(product_part)

            if product_mol and all(r for r in reactants_mols):
                # Look for hydroxyl in reactants
                oh_pattern = Chem.MolFromSmarts("[OH]")
                benzyl_pattern = Chem.MolFromSmarts("[CH2][c]1[cH][cH][cH][cH][cH]1")

                # Look for ether in product
                ether_pattern = Chem.MolFromSmarts("[O][CH2][c]")

                reactants_have_oh = any(r.HasSubstructMatch(oh_pattern) for r in reactants_mols)
                reactants_have_benzyl = any(
                    r.HasSubstructMatch(benzyl_pattern) for r in reactants_mols
                )
                product_has_ether = product_mol.HasSubstructMatch(ether_pattern)

                if reactants_have_oh and reactants_have_benzyl and product_has_ether:
                    print(f"O-alkylation detected at depth {depth}")
                    o_alkylation_count += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return o_alkylation_count >= 2  # Return True if at least 2 O-alkylations are detected
