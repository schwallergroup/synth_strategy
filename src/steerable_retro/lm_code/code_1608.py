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
    This function detects if the synthetic route involves azide chemistry,
    specifically looking for azide introduction and transformation.
    """
    azide_found = False

    def dfs_traverse(node):
        nonlocal azide_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Check for azide pattern in reactants or products
                azide_pattern = Chem.MolFromSmarts("[N-]=[N+]=N")

                try:
                    reactant_mol = Chem.MolFromSmiles(reactants)
                    product_mol = Chem.MolFromSmiles(product)

                    if reactant_mol and product_mol:
                        if reactant_mol.HasSubstructMatch(
                            azide_pattern
                        ) or product_mol.HasSubstructMatch(azide_pattern):
                            print(f"Found azide pattern in reaction: {rsmi}")
                            azide_found = True
                except:
                    print(f"Error processing SMILES: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return azide_found
