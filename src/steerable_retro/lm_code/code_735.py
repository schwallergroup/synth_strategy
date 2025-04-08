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
    Detects synthesis strategy involving introduction of allyloxy group.
    """
    allyloxy_introduction = False

    def dfs_traverse(node):
        nonlocal allyloxy_introduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for allyloxy group introduction
            allyl_pattern = Chem.MolFromSmarts("[CH2]=[CH]-[CH2]-[Br,Cl,I]")
            phenol_pattern = Chem.MolFromSmarts("[c]-[OH]")
            allyloxy_pattern = Chem.MolFromSmarts("[CH2]=[CH]-[CH2]-[O]-[c]")

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(allyloxy_pattern):
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and (
                        reactant_mol.HasSubstructMatch(allyl_pattern)
                        or reactant_mol.HasSubstructMatch(phenol_pattern)
                    ):
                        print("Detected allyloxy group introduction")
                        allyloxy_introduction = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return allyloxy_introduction
