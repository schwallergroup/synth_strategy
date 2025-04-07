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
    This function detects if the synthesis involves formation of a 5-membered
    heterocyclic ring (specifically imidazolidinedione/hydantoin).
    """
    hydantoin_pattern = Chem.MolFromSmarts("[#6]1[#7][#6](=[#8])[#7][#6](=[#8])1")
    found_ring_formation = False

    def dfs_traverse(node):
        nonlocal found_ring_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            # Check if product contains the heterocycle but reactants don't
            if product_mol and product_mol.HasSubstructMatch(hydantoin_pattern):
                if not any(
                    mol and mol.HasSubstructMatch(hydantoin_pattern) for mol in reactant_mols
                ):
                    print("Found heterocycle formation reaction")
                    found_ring_formation = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_ring_formation
