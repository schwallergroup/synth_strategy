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
    Detects if the synthetic route contains an amide formation step
    (carboxylic acid + amine -> amide)
    """
    found_amide_formation = False

    def dfs_traverse(node):
        nonlocal found_amide_formation

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check for carboxylic acid in reactants
            acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")

            # Check for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[N;H]")

            # Check for amide in product
            amide_pattern = Chem.MolFromSmarts("C(=O)[N]")

            reactants = reactants_part.split(".")
            has_acid = False
            has_amine = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(acid_pattern):
                        has_acid = True
                    if mol and mol.HasSubstructMatch(amine_pattern):
                        has_amine = True
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product_part)
                has_amide = product_mol and product_mol.HasSubstructMatch(amide_pattern)
            except:
                has_amide = False

            if has_acid and has_amine and has_amide:
                print("Found amide formation: carboxylic acid + amine -> amide")
                found_amide_formation = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_amide_formation
