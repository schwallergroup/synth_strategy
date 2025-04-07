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
    Detects aldehyde protection/deprotection sequence in the route
    """
    acetal_found = False
    acetal_to_aldehyde = False

    def dfs_traverse(node):
        nonlocal acetal_found, acetal_to_aldehyde

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for acetal in reactants
            acetal_pattern = Chem.MolFromSmarts("[#6]([#8][#6])([#8][#6])[#6]")
            aldehyde_pattern = Chem.MolFromSmarts("[#6](=[#8])[#1,#6]")

            # Check for acetal in reactants
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(acetal_pattern):
                        acetal_found = True
                        print(f"Found acetal in reactant: {reactant}")

                        # Check if product has aldehyde
                        try:
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol and product_mol.HasSubstructMatch(
                                aldehyde_pattern
                            ):
                                acetal_to_aldehyde = True
                                print(f"Found acetal to aldehyde transformation")
                        except:
                            pass
                except:
                    continue

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    result = acetal_found and acetal_to_aldehyde
    print(f"Aldehyde protection/deprotection detected: {result}")
    return result
