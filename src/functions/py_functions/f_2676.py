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
    Detects if the synthesis route involves sequential modifications of aromatic rings,
    specifically looking for hydroxyl to chloro and hydroxyl to methoxy transformations.
    """
    hydroxyl_to_chloro_found = False
    hydroxyl_to_methoxy_found = False

    def dfs_traverse(node):
        nonlocal hydroxyl_to_chloro_found, hydroxyl_to_methoxy_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for hydroxyl to chloro conversion
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    product_mol = Chem.MolFromSmiles(product)

                    if not (reactant_mol and product_mol):
                        continue

                    # Check for hydroxyl in reactant
                    hydroxyl_pattern = Chem.MolFromSmarts("[O][c]")
                    if reactant_mol.HasSubstructMatch(hydroxyl_pattern):

                        # Check for chloro in product
                        chloro_pattern = Chem.MolFromSmarts("[Cl][c]")
                        if product_mol.HasSubstructMatch(chloro_pattern):
                            hydroxyl_to_chloro_found = True
                            print("Hydroxyl to chloro conversion detected")

                        # Check for methoxy in product
                        methoxy_pattern = Chem.MolFromSmarts("[C][O][c]")
                        if product_mol.HasSubstructMatch(methoxy_pattern):
                            hydroxyl_to_methoxy_found = True
                            print("Hydroxyl to methoxy conversion detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    # Return True if both modifications are found
    return hydroxyl_to_chloro_found and hydroxyl_to_methoxy_found
