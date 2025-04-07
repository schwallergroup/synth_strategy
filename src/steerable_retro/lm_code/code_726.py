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
    This function detects a synthetic strategy involving piperazine ring formation
    while maintaining Boc protection throughout the synthesis.
    """
    piperazine_formed = False
    boc_protected = False

    def dfs_traverse(node):
        nonlocal piperazine_formed, boc_protected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for piperazine formation
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    piperazine_pattern = Chem.MolFromSmarts("[N]1[C][C][N][C][C]1")
                    if product_mol.HasSubstructMatch(piperazine_pattern):
                        # Check if any reactant doesn't have piperazine
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and not reactant_mol.HasSubstructMatch(
                                piperazine_pattern
                            ):
                                piperazine_formed = True
                                print("Detected piperazine ring formation")
                                break

                # Check for Boc protection
                boc_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]([C])([C])[C]")
                if product_mol and product_mol.HasSubstructMatch(boc_pattern):
                    boc_protected = True
                    print("Detected Boc protection")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return piperazine_formed and boc_protected
