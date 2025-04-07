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
    This function detects if the synthesis includes a morpholine amide coupling step.
    """
    has_morpholine_coupling = False

    def dfs_traverse(node):
        nonlocal has_morpholine_coupling

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for morpholine pattern in reactants
            morpholine_pattern = Chem.MolFromSmarts("[N]1[CH2][CH2][O][CH2][CH2]1")
            amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

            product_mol = Chem.MolFromSmiles(product)

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(morpholine_pattern):
                    # Check if product has amide bond
                    if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                        has_morpholine_coupling = True
                        print("Detected morpholine amide coupling")
                        break

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_morpholine_coupling
