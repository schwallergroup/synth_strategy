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
    This function detects late-stage alcohol activation via mesylation.
    """
    found_mesylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_mesylation

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Check only late-stage reactions (depth 0 or 1)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants.split(".")]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    # Define SMARTS patterns
                    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
                    mesylate_pattern = Chem.MolFromSmarts(
                        "[SX4](=[OX1])(=[OX1])([OX2])[#6]"
                    )

                    # Check for mesylation (alcohol → mesylate)
                    alcohol_in_reactants = any(
                        r and r.HasSubstructMatch(alcohol_pattern)
                        for r in reactants_mols
                        if r
                    )
                    mesylate_in_product = product_mol.HasSubstructMatch(
                        mesylate_pattern
                    )

                    if alcohol_in_reactants and mesylate_in_product:
                        print(f"Detected late-stage mesylation at depth {depth}")
                        found_mesylation = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_mesylation
