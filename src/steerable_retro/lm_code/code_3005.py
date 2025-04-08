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
    This function detects a halogen exchange strategy (Cl to I) on a heterocycle.
    """
    halogen_exchange_detected = False

    def dfs_traverse(node):
        nonlocal halogen_exchange_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for chloro-heterocycle to iodo-heterocycle transformation
                chloro_heterocycle_pattern = Chem.MolFromSmarts("c1[n]c(Cl)[n]cc1")
                iodo_heterocycle_pattern = Chem.MolFromSmarts("c1[n]c(I)[n]cc1")

                try:
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            chloro_heterocycle_pattern
                        ):
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol and product_mol.HasSubstructMatch(
                                iodo_heterocycle_pattern
                            ):
                                halogen_exchange_detected = True
                                print(
                                    "Halogen exchange detected: Chloro-heterocycle -> Iodo-heterocycle"
                                )
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return halogen_exchange_detected
