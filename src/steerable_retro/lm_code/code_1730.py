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
    This function detects if the synthesis route involves formation of an ether by replacing
    a halogen (particularly fluorine) with an alkoxy group.
    """
    halogen_to_ether = False

    def dfs_traverse(node):
        nonlocal halogen_to_ether

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactant contains aryl fluoride and product contains aryl ether
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("cF")):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("cOC")):
                            halogen_to_ether = True
                            print(f"Halogen to ether transformation detected: {rsmi}")
                            break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return halogen_to_ether
