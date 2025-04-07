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
    This function detects if the synthetic route involves transformation
    of a chloride to an azide group.
    """
    transformation_found = False

    def dfs_traverse(node):
        nonlocal transformation_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Check for chloride in reactants
                chloride_pattern = Chem.MolFromSmarts("[#6]-[Cl]")
                # Check for azide in products
                azide_pattern = Chem.MolFromSmarts("[#6]-[N]=[N+]=[N-]")

                try:
                    reactant_mol = Chem.MolFromSmiles(reactants)
                    product_mol = Chem.MolFromSmiles(product)

                    if reactant_mol and product_mol:
                        if reactant_mol.HasSubstructMatch(
                            chloride_pattern
                        ) and product_mol.HasSubstructMatch(azide_pattern):
                            print(
                                f"Found chloride to azide transformation in reaction: {rsmi}"
                            )
                            transformation_found = True
                except:
                    print(f"Error processing SMILES: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return transformation_found
