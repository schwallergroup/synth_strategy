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
    Detects if the synthesis route uses etherification via tosylate displacement.
    """
    etherification_detected = False
    tosylate_pattern = Chem.MolFromSmarts("[#6]c1ccc(cc1)[S](=[#8])(=[#8])[#8][#6]")
    ether_pattern = Chem.MolFromSmarts("[#6][#8][#6]")

    def dfs_traverse(node):
        nonlocal etherification_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for tosylate in reactants and ether in product
                tosylate_in_reactants = any(
                    Chem.MolFromSmiles(r).HasSubstructMatch(tosylate_pattern)
                    if Chem.MolFromSmiles(r)
                    else False
                    for r in reactants
                )
                phenol_in_reactants = any("[OH]c" in r for r in reactants)

                product_mol = Chem.MolFromSmiles(product)
                if (
                    product_mol
                    and product_mol.HasSubstructMatch(ether_pattern)
                    and (tosylate_in_reactants or phenol_in_reactants)
                ):
                    print(
                        f"Found etherification at depth: {node.get('metadata', {}).get('depth', 'unknown')}"
                    )
                    etherification_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Etherification via tosylate strategy detected: {etherification_detected}")
    return etherification_detected
