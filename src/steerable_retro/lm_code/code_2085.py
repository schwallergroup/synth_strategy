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
    Detects if the route contains a transformation from alcohol to mesylate.
    """
    mesylation_found = False

    def dfs_traverse(node):
        nonlocal mesylation_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol in reactants
            has_alcohol = False
            for reactant in reactants:
                if reactant:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        alcohol_pattern = Chem.MolFromSmarts("[C][OH]")
                        if mol.HasSubstructMatch(alcohol_pattern):
                            has_alcohol = True

            # Check for mesylate in product
            if has_alcohol:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol:
                    mesylate_pattern = Chem.MolFromSmarts("[C][O][S](=O)(=O)[C]")
                    if prod_mol.HasSubstructMatch(mesylate_pattern):
                        print("Found alcohol to mesylate transformation")
                        mesylation_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return mesylation_found
