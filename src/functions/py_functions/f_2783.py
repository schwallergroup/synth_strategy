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
    Detects a strategy involving the formation of multiple heterocycles
    (specifically thiazole and isoxazoline rings).
    """
    thiazole_formation = False
    isoxazoline_formation = False

    def dfs_traverse(node):
        nonlocal thiazole_formation, isoxazoline_formation

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol:
                # Check for thiazole formation
                thiazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#6][#16][#6]1")
                if product_mol.HasSubstructMatch(thiazole_pattern):
                    # Check if thiazole wasn't in reactants
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    if not any(
                        mol and mol.HasSubstructMatch(thiazole_pattern)
                        for mol in reactant_mols
                    ):
                        print("Detected thiazole ring formation")
                        thiazole_formation = True

                # Check for isoxazoline formation
                isoxazoline_pattern = Chem.MolFromSmarts("[#6]1[#6][#8][#7]=[#6]1")
                if product_mol.HasSubstructMatch(isoxazoline_pattern):
                    # Check if isoxazoline wasn't in reactants
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    if not any(
                        mol and mol.HasSubstructMatch(isoxazoline_pattern)
                        for mol in reactant_mols
                    ):
                        print("Detected isoxazoline ring formation")
                        isoxazoline_formation = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thiazole_formation and isoxazoline_formation
