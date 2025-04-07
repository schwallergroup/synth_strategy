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
    Detects if the synthetic route involves a convergent step where two complex fragments
    are joined via ester formation.
    """
    convergent_ester_found = False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_ester_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an ester formation reaction
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                ester_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#8]-[#6]")
                if product_mol.HasSubstructMatch(ester_pattern):
                    # Check if we have at least two complex reactants
                    complex_reactants = 0
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if (
                            reactant_mol and reactant_mol.GetNumAtoms() > 6
                        ):  # Arbitrary threshold for "complex"
                            complex_reactants += 1

                    if complex_reactants >= 2:
                        convergent_ester_found = True
                        print(f"Convergent ester formation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return convergent_ester_found
