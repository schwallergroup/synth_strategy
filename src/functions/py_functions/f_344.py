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
    This function detects a strategy involving multiple alcohol to leaving group
    conversions (mesylate, chloride, etc.) throughout the synthesis.
    """
    alcohol_conversions = 0

    def dfs_traverse(node):
        nonlocal alcohol_conversions

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol in reactants
            alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")

            # Check for mesylate or chloride in product
            mesylate_pattern = Chem.MolFromSmarts("[#6]-[#8]-S(=O)(=O)-[#6]")
            chloride_pattern = Chem.MolFromSmarts("[#6]-[Cl]")

            try:
                product_mol = Chem.MolFromSmiles(product)

                for reactant in reactants:
                    if reactant:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            alcohol_pattern
                        ):
                            if product_mol and (
                                product_mol.HasSubstructMatch(mesylate_pattern)
                                or product_mol.HasSubstructMatch(chloride_pattern)
                            ):
                                print(
                                    "Detected alcohol to leaving group conversion in reaction:",
                                    rsmi,
                                )
                                alcohol_conversions += 1
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return alcohol_conversions >= 2  # At least 2 alcohol conversions
