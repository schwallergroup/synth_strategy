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
    Detects if the route incorporates a cyclopropyl group.
    """
    found_cyclopropyl = False

    def dfs_traverse(node):
        nonlocal found_cyclopropyl

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Pattern for cyclopropyl group
            cyclopropyl_pattern = Chem.MolFromSmarts("[C]1[C][C]1")

            try:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(cyclopropyl_pattern):
                    # Check if any reactant has cyclopropyl
                    has_cyclopropyl_reactant = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            cyclopropyl_pattern
                        ):
                            has_cyclopropyl_reactant = True
                            break

                    if has_cyclopropyl_reactant:
                        found_cyclopropyl = True
                        print(
                            f"Found cyclopropyl incorporation at depth {node.get('depth', 'unknown')}"
                        )
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_cyclopropyl
