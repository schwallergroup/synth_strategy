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
    Detects the formation of an amide bond with a cyclopropylamine.
    """
    found_cyclopropyl_amide = False

    def dfs_traverse(node):
        nonlocal found_cyclopropyl_amide

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Patterns
            cyclopropyl_amine_pattern = Chem.MolFromSmarts("[NH2][C]1[C][C]1")
            carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
            cyclopropyl_amide_pattern = Chem.MolFromSmarts("[C](=[O])[NH][C]1[C][C]1")

            # Check if reactants include cyclopropylamine and carboxylic acid
            has_cyclopropyl_amine = False
            has_carboxylic_acid = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(cyclopropyl_amine_pattern):
                        has_cyclopropyl_amine = True
                    if mol and mol.HasSubstructMatch(carboxylic_acid_pattern):
                        has_carboxylic_acid = True
                except:
                    continue

            # Check if product has cyclopropyl amide
            try:
                product_mol = Chem.MolFromSmiles(product)
                if (
                    product_mol
                    and product_mol.HasSubstructMatch(cyclopropyl_amide_pattern)
                    and has_cyclopropyl_amine
                    and has_carboxylic_acid
                ):
                    print("Found cyclopropyl amide formation")
                    found_cyclopropyl_amide = True
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_cyclopropyl_amide
