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
    This function detects if the synthesis route involves transformation of a halogen (particularly iodine)
    to a nitrile group.
    """
    halogen_to_nitrile = False

    def dfs_traverse(node):
        nonlocal halogen_to_nitrile

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactant contains iodine and product contains nitrile
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("cI")):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("cC#N")
                        ):
                            halogen_to_nitrile = True
                            print(f"Halogen to nitrile transformation detected: {rsmi}")
                            break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return halogen_to_nitrile
