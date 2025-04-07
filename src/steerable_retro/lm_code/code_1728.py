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
    This function detects if the synthesis route involves sequential functionalization
    of a pyridine core scaffold.
    """
    pyridine_modifications = 0

    def dfs_traverse(node):
        nonlocal pyridine_modifications

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains pyridine
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("c1ccncc1")):
                    # Check if this is a modification of the pyridine ring
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("c1ccncc1")
                        ):
                            pyridine_modifications += 1
                            print(f"Pyridine modification detected: {rsmi}")
                            break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return pyridine_modifications >= 2  # At least 2 modifications to consider it a strategy
