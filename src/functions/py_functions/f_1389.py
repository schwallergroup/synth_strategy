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
    Detects if the synthesis route involves formation of a thiocarbamate group (O-C(=S)-O)
    in the late stage (low depth) of the synthesis.
    """
    thiocarbamate_formed = False
    formation_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal thiocarbamate_formed, formation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if thiocarbamate is formed in this reaction
            product_mol = Chem.MolFromSmiles(product_smiles)
            thiocarbamate_pattern = Chem.MolFromSmarts("[O]-[C](=[S])-[O]")

            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and not reactant_mol.HasSubstructMatch(
                    thiocarbamate_pattern
                ):
                    if product_mol and product_mol.HasSubstructMatch(
                        thiocarbamate_pattern
                    ):
                        thiocarbamate_formed = True
                        formation_depth = min(formation_depth, depth)
                        print(f"Thiocarbamate formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if thiocarbamate was formed in late stage (depth <= 1)
    late_stage = formation_depth <= 1

    if thiocarbamate_formed and late_stage:
        print(f"Late-stage thiocarbamate formation detected at depth {formation_depth}")
        return True
    return False
