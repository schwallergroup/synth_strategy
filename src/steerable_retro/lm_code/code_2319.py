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
    Detects if the synthesis route incorporates morpholine rings via nucleophilic substitution.
    """
    morpholine_addition = False

    def dfs_traverse(node, depth=0):
        nonlocal morpholine_addition

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Morpholine pattern
                morpholine_pattern = Chem.MolFromSmarts("[#8]1[#6][#6][#7][#6][#6]1")

                # Check if morpholine is in reactants
                morpholine_in_reactants = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(morpholine_pattern):
                        morpholine_in_reactants = True
                        break

                # Check if morpholine is in product
                prod_mol = Chem.MolFromSmiles(product)
                morpholine_in_product = prod_mol and prod_mol.HasSubstructMatch(morpholine_pattern)

                # If morpholine is in both reactants and product, and we have a halide in reactants,
                # it's likely a nucleophilic substitution
                if morpholine_in_reactants and morpholine_in_product:
                    halide_pattern = Chem.MolFromSmarts("[#6]-[#17,#35,#53]")
                    halide_in_reactants = False

                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(halide_pattern):
                            halide_in_reactants = True
                            break

                    if halide_in_reactants:
                        morpholine_addition = True
                        print(
                            "Detected morpholine incorporation via nucleophilic substitution at depth",
                            depth,
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return morpholine_addition
