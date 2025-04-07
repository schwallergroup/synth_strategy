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
    Detects the transformation of a nitrile group to an amidoxime group.
    """
    transformation_detected = False

    def dfs_traverse(node):
        nonlocal transformation_detected

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitrile in reactants
            nitrile_pattern = Chem.MolFromSmarts("[#6]#[#7]")
            nitrile_present = False
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(nitrile_pattern):
                        nitrile_present = True
                        print(f"Found nitrile in reactant: {reactant}")
                except:
                    continue

            # Check for amidoxime in product
            amidoxime_pattern = Chem.MolFromSmarts("[#6](=[#7])[#7][#8]")
            amidoxime_present = False
            try:
                mol = Chem.MolFromSmiles(product)
                if mol and mol.HasSubstructMatch(amidoxime_pattern):
                    amidoxime_present = True
                    print(f"Found amidoxime in product: {product}")
            except:
                pass

            # If both conditions are met, we've found our strategy
            if nitrile_present and amidoxime_present:
                transformation_detected = True
                print("Detected nitrile to amidoxime transformation")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return transformation_detected
