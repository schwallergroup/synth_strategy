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
    Detects if the synthesis incorporates a trifluoromethyl-containing aromatic group.
    """
    cf3_incorporation = False

    def dfs_traverse(node):
        nonlocal cf3_incorporation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if CF3 group is present in reactants but not in product
            cf3_pattern = Chem.MolFromSmarts("[#6]-[#6]([F])([F])[F]")

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                product_mol = Chem.MolFromSmiles(product_part)

                if product_mol and any(
                    r and r.HasSubstructMatch(cf3_pattern) for r in reactants if r
                ):
                    if product_mol.HasSubstructMatch(cf3_pattern):
                        print("Detected incorporation of trifluoromethyl-containing aromatic group")
                        cf3_incorporation = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return cf3_incorporation
