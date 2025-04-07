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
    This function detects if the route contains aromatization of a saturated heterocycle,
    specifically tetrahydropyridine to pyridine conversion.
    """
    found = False

    def dfs_traverse(node):
        nonlocal found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for tetrahydropyridine pattern in reactants
            reactant_mol = Chem.MolFromSmiles(reactants)
            tetrahydropyridine_pattern = Chem.MolFromSmarts(
                "[#6]1[#6][#6][#7][#6][#6]1"
            )

            # Check for pyridine pattern in product
            product_mol = Chem.MolFromSmiles(product)
            pyridine_pattern = Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#7]:[#6]:[#6]:1")

            if reactant_mol and product_mol:
                if reactant_mol.HasSubstructMatch(
                    tetrahydropyridine_pattern
                ) and product_mol.HasSubstructMatch(pyridine_pattern):
                    print(
                        "Found heterocycle aromatization (tetrahydropyridine to pyridine)"
                    )
                    found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return found
