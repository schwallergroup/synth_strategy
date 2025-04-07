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
    Detects if the synthesis involves construction of a morpholine core
    through a ring-forming reaction.
    """
    found_morpholine_formation = False

    def dfs_traverse(node):
        nonlocal found_morpholine_formation

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            # Look for morpholine ring in product
            morpholine_pattern = Chem.MolFromSmarts(
                "[NX3]1-[CX4]-[CX4]-[OX2]-[CX4]-[CX4]-1"
            )

            # Check if morpholine is in product but not in all reactants
            has_morpholine_product = (
                product_mol is not None
                and product_mol.HasSubstructMatch(morpholine_pattern)
            )
            all_reactants_have_morpholine = all(
                mol is not None and mol.HasSubstructMatch(morpholine_pattern)
                for mol in reactant_mols
                if mol is not None
            )

            if has_morpholine_product and not all_reactants_have_morpholine:
                found_morpholine_formation = True
                print("Found morpholine ring formation step")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_morpholine_formation
