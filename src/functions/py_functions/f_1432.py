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
    This function detects a synthetic strategy involving Boc protection of an amine.
    """
    has_boc_protection = False

    def dfs_traverse(node):
        nonlocal has_boc_protection

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc protection
                boc_pattern = Chem.MolFromSmarts("[#6]OC(=O)N")
                boc_reagent_pattern = "CC(C)(C)OC(=O)O"

                if any(boc_reagent_pattern in r for r in reactants):
                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(boc_pattern):
                            has_boc_protection = True
                            print("Found Boc protection")
                    except:
                        pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return has_boc_protection
