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
    Detects if the synthesis route involves transformation of a nitrile group to an amidine group.
    """
    nitrile_transformed = False

    def dfs_traverse(node):
        nonlocal nitrile_transformed

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check for nitrile in reactants
                    nitrile_pattern = Chem.MolFromSmarts("[#6]#[#7]")
                    amidine_pattern = Chem.MolFromSmarts("[#6](=[#7])[#7]")

                    reactant_has_nitrile = False
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(nitrile_pattern):
                            reactant_has_nitrile = True
                            break

                    # Check for amidine in product
                    product_mol = Chem.MolFromSmiles(product)
                    product_has_amidine = product_mol and product_mol.HasSubstructMatch(
                        amidine_pattern
                    )

                    if reactant_has_nitrile and product_has_amidine:
                        nitrile_transformed = True
                        print(f"Nitrile to amidine transformation detected in reaction: {rsmi}")
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitrile to amidine transformation detected: {nitrile_transformed}")
    return nitrile_transformed
