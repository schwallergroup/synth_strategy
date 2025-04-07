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
    This function detects functionalization at the C2 position of indole.
    """
    has_c2_functionalization = False

    def dfs_traverse(node):
        nonlocal has_c2_functionalization

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for C2 functionalization of indole
                indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")
                c2_functionalized_indole = Chem.MolFromSmarts("c1ccc2[nH]c(C(=O)[OH])cc2c1")

                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(c2_functionalized_indole):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(indole_pattern):
                            if not reactant_mol.HasSubstructMatch(c2_functionalized_indole):
                                has_c2_functionalization = True
                                print("Detected indole C2 functionalization")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_c2_functionalization
