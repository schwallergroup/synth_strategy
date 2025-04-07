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
    This function detects the formation of a guanidine-like moiety (NH-C(=NH)-Ar).
    """
    guanidine_formation_found = False

    def dfs_traverse(node):
        nonlocal guanidine_formation_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for guanidine pattern in product but not in reactants
            guanidine_pattern = "[#7]-C(=[#7])-[#6]"

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts(guanidine_pattern)
            ):
                # Check if reactants don't have this pattern
                reactants_have_pattern = False
                for r in reactants:
                    r_mol = Chem.MolFromSmiles(r)
                    if r_mol and r_mol.HasSubstructMatch(
                        Chem.MolFromSmarts(guanidine_pattern)
                    ):
                        reactants_have_pattern = True
                        break

                if not reactants_have_pattern:
                    guanidine_formation_found = True
                    print("Guanidine-like moiety formation detected")

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return guanidine_formation_found
