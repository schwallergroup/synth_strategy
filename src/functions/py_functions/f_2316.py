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
    This function detects if the synthesis includes sequential deprotection of different protecting groups.
    """
    deprotection_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal deprotection_depths

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc deprotection
            boc_pattern = Chem.MolFromSmarts(
                "[#6]([#6])([#6])([#6])[#8][#6](=[#8])[#7]"
            )
            has_boc_in_reactants = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(boc_pattern):
                    has_boc_in_reactants = True
                    break

            product_mol = Chem.MolFromSmiles(product)
            has_boc_in_product = product_mol and product_mol.HasSubstructMatch(
                boc_pattern
            )

            if has_boc_in_reactants and not has_boc_in_product:
                deprotection_depths.append((depth, "Boc"))

            # Check for methoxy deprotection
            methoxy_pattern = Chem.MolFromSmarts("[#6][#8][#6]:[#6]")
            hydroxy_pattern = Chem.MolFromSmarts("[#8][#6]:[#6]")

            has_methoxy_in_reactants = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(methoxy_pattern):
                    has_methoxy_in_reactants = True
                    break

            has_hydroxy_in_product = product_mol and product_mol.HasSubstructMatch(
                hydroxy_pattern
            )
            has_methoxy_in_product = product_mol and product_mol.HasSubstructMatch(
                methoxy_pattern
            )

            if (
                has_methoxy_in_reactants
                and has_hydroxy_in_product
                and not has_methoxy_in_product
            ):
                deprotection_depths.append((depth, "Methoxy"))

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Sort deprotection events by depth
    deprotection_depths.sort()

    # Check if we have at least 2 different deprotection events
    if len(deprotection_depths) >= 2:
        # Check if they are different types
        deprotection_types = set([d[1] for d in deprotection_depths])
        if len(deprotection_types) >= 2:
            print(f"Detected sequential deprotection strategy: {deprotection_depths}")
            return True

    return False
