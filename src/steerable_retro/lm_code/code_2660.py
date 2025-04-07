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
    This function detects if the synthesis involves functionalization of a nitrogen heterocycle core.
    """
    has_nitrogen_heterocycle = False
    functionalization_steps = 0

    def dfs_traverse(node):
        nonlocal has_nitrogen_heterocycle, functionalization_steps

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)

                # Check if product contains nitrogen heterocycle
                # Common nitrogen heterocycles: pyridine, quinoline, indole, etc.
                nitrogen_heterocycle_patterns = [
                    Chem.MolFromSmarts("c1cccnc1"),  # pyridine
                    Chem.MolFromSmarts("c1ccc2ncccc2c1"),  # quinoline
                    Chem.MolFromSmarts("c1ccc2[nH]ccc2c1"),  # indole
                    Chem.MolFromSmarts("c1ccc2c(c1)cccn2"),  # alternative quinoline
                ]

                if product_mol:
                    for pattern in nitrogen_heterocycle_patterns:
                        if product_mol.HasSubstructMatch(pattern):
                            has_nitrogen_heterocycle = True

                            # Check if this is a functionalization step
                            for reactant in reactants:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol and reactant_mol.HasSubstructMatch(pattern):
                                    # If both reactant and product have the heterocycle, and product has more atoms,
                                    # it's likely a functionalization
                                    if product_mol.GetNumAtoms() > reactant_mol.GetNumAtoms():
                                        functionalization_steps += 1
                                        print(
                                            f"Detected nitrogen heterocycle functionalization: {rsmi}"
                                        )
                                    break
            except:
                print(f"Error in processing molecules for nitrogen heterocycle check: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Return True if there's a nitrogen heterocycle with at least 2 functionalization steps
    result = has_nitrogen_heterocycle and functionalization_steps >= 2
    print(f"Nitrogen heterocycle core functionalization: {result}")
    print(f"  - Has nitrogen heterocycle: {has_nitrogen_heterocycle}")
    print(f"  - Functionalization steps: {functionalization_steps}")

    return result
