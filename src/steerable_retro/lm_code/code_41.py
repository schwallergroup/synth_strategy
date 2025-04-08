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
    This function detects if the synthesis route involves N-Boc protection strategy.
    """
    boc_protection_detected = False

    def dfs_traverse(node):
        nonlocal boc_protection_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for primary/secondary amine in reactants
            amine_pattern = Chem.MolFromSmarts("[N;H1,H2]")
            # Check for Boc-protected amine in product
            boc_pattern = Chem.MolFromSmarts("[N][C](=[O])[O][C]([C])([C])[C]")

            for reactant in reactants:
                if reactant and Chem.MolFromSmiles(reactant) and amine_pattern:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(amine_pattern):
                        if product and Chem.MolFromSmiles(product) and boc_pattern:
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol and product_mol.HasSubstructMatch(boc_pattern):
                                boc_protection_detected = True
                                print("Detected N-Boc protection")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return boc_protection_detected
