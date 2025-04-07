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
    This function detects a strategy involving thiazole ring formation from a thiourea intermediate.
    """
    thiazole_formed = False
    thiourea_present = False

    def dfs_traverse(node):
        nonlocal thiazole_formed, thiourea_present

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains thiazole
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    thiazole_pattern = Chem.MolFromSmarts("c1nc([#6])cs1")
                    if product_mol.HasSubstructMatch(thiazole_pattern):
                        # Check if reactants contain thiourea
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                thiourea_pattern = Chem.MolFromSmarts(
                                    "[#6]-[NH]-[C](=[S])-[NH2]"
                                )
                                if reactant_mol.HasSubstructMatch(thiourea_pattern):
                                    thiourea_present = True
                                    thiazole_formed = True
                                    print("Found thiazole formation from thiourea")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thiazole_formed and thiourea_present
