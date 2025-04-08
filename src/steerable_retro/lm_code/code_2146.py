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
    This function detects the strategy of forming a tetrazole ring from a nitrile group.
    """
    tetrazole_formed = False

    def dfs_traverse(node):
        nonlocal tetrazole_formed

        if node["type"] == "reaction":
            # Check if this is a reaction node with metadata
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitrile in reactants
                nitrile_pattern = Chem.MolFromSmarts("[#6]#[#7]")
                for reactant in reactants:
                    if reactant:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(nitrile_pattern):
                                # Check for tetrazole in product
                                tetrazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#7][#7][#7]1")
                                product_mol = Chem.MolFromSmiles(product)
                                if product_mol and product_mol.HasSubstructMatch(tetrazole_pattern):
                                    print("Detected tetrazole formation from nitrile")
                                    tetrazole_formed = True
                        except:
                            continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return tetrazole_formed
