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
    Detects if the synthesis involves an amide bond disconnection.
    """
    amide_disconnection_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_disconnection_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product)

            if product_mol:
                # Check for amide pattern in product
                amide_pattern = Chem.MolFromSmarts("[NH][C](=[O])")
                if product_mol.HasSubstructMatch(amide_pattern):
                    # Check if reactants contain acid/acid derivative and amine
                    acid_pattern = Chem.MolFromSmarts("[C](=[O])[O,Cl,Br,I]")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")

                    has_acid = any(
                        Chem.MolFromSmiles(r)
                        and Chem.MolFromSmiles(r).HasSubstructMatch(acid_pattern)
                        for r in reactants
                        if r
                    )
                    has_amine = any(
                        Chem.MolFromSmiles(r)
                        and Chem.MolFromSmiles(r).HasSubstructMatch(amine_pattern)
                        for r in reactants
                        if r
                    )

                    if has_acid and has_amine:
                        amide_disconnection_detected = True
                        print(f"Amide disconnection detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return amide_disconnection_detected
