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
    This function detects a strategy involving MEM protection and deprotection of a sulfonamide nitrogen.
    """
    has_mem_protection = False
    has_mem_deprotection = False

    def dfs_traverse(node):
        nonlocal has_mem_protection, has_mem_deprotection

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for MEM protection (N-H to N-MEM)
            if (
                any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(
                        Chem.MolFromSmarts("[NH][S](=[O])(=[O])[c]")
                    )
                    for r in reactants
                )
                and Chem.MolFromSmiles(product)
                and Chem.MolFromSmiles(product).HasSubstructMatch(
                    Chem.MolFromSmarts(
                        "[N]([CH2][O][CH2][CH2][O][CH3])[S](=[O])(=[O])[c]"
                    )
                )
            ):
                has_mem_protection = True
                print("Found MEM protection reaction")

            # Check for MEM deprotection (N-MEM to N-H)
            if any(
                Chem.MolFromSmiles(r)
                and Chem.MolFromSmiles(r).HasSubstructMatch(
                    Chem.MolFromSmarts(
                        "[N]([CH2][O][CH2][CH2][O][CH3])[S](=[O])(=[O])[c]"
                    )
                )
                for r in reactants
            ) and any(
                Chem.MolFromSmiles(p)
                and Chem.MolFromSmiles(p).HasSubstructMatch(
                    Chem.MolFromSmarts("[NH][S](=[O])(=[O])[c]")
                )
                for p in product.split(".")
            ):
                has_mem_deprotection = True
                print("Found MEM deprotection reaction")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return has_mem_protection and has_mem_deprotection
