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
    Detects Fischer-type indole synthesis from nitro aromatic precursor.
    Looks for a reaction that converts a nitro-containing aromatic to an indole structure.
    """
    indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")
    nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])-[O-]")

    found_indole_formation = False

    def dfs_traverse(node):
        nonlocal found_indole_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains indole
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(indole_pattern):
                    # Check if any reactant contains nitro group
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            nitro_pattern
                        ):
                            print("Found indole formation from nitro compound")
                            found_indole_formation = True
                            break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_indole_formation
