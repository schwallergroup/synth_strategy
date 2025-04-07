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
    This function detects a synthetic strategy involving formation of a heterocyclic ring
    via nitro reduction and cyclization.
    """
    ring_formation_found = False

    def dfs_traverse(node):
        nonlocal ring_formation_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for benzoxazinone formation
            benzoxazinone_pattern = Chem.MolFromSmarts(
                "[#8]1[#6][#6](=O)[#7][#6]2[#6][#6][#6][#6][#6]12"
            )
            product_mol = Chem.MolFromSmiles(product)

            if (
                product_mol
                and benzoxazinone_pattern
                and product_mol.HasSubstructMatch(benzoxazinone_pattern)
            ):
                # Check if reactants contain nitro group
                nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if (
                        reactant_mol
                        and nitro_pattern
                        and reactant_mol.HasSubstructMatch(nitro_pattern)
                    ):
                        print(
                            "Found heterocyclic ring formation via nitro reduction and cyclization"
                        )
                        ring_formation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return ring_formation_found
