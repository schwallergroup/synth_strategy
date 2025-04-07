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
    This function detects if the synthesis route involves transformation of an alcohol to an ester.
    """
    transformation_found = False

    def dfs_traverse(node):
        nonlocal transformation_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol to ester transformation
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                ester_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]=[#8]")
                if product_mol.HasSubstructMatch(ester_pattern):
                    # Check if any reactant has an alcohol
                    alcohol_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")
                            if reactant_mol.HasSubstructMatch(alcohol_pattern):
                                alcohol_in_reactants = True

                    # If ester is in product and alcohol in reactants, it's an alcohol to ester transformation
                    if alcohol_in_reactants:
                        transformation_found = True
                        print("Alcohol to ester transformation detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return transformation_found
