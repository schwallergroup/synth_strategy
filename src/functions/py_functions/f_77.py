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
    This function detects O-deprotection reactions, particularly methoxy to hydroxyl conversion.
    """
    o_deprotection_detected = False

    def dfs_traverse(node):
        nonlocal o_deprotection_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Patterns for O-deprotection
                methoxy_pattern = Chem.MolFromSmarts("[#6][O][CH3]")  # Methoxy group
                hydroxyl_pattern = Chem.MolFromSmarts("[#6][OH]")  # Hydroxyl group

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and product_mol.HasSubstructMatch(hydroxyl_pattern):
                    if any(
                        r and r.HasSubstructMatch(methoxy_pattern)
                        for r in reactant_mols
                    ):
                        print("Detected O-deprotection: methoxy to hydroxyl conversion")
                        o_deprotection_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return o_deprotection_detected
