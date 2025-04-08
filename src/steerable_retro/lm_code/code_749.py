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
    Detects a synthetic strategy involving tosylation of an alcohol (leaving group installation).
    """
    has_alcohol_tosylation = False

    def dfs_traverse(node):
        nonlocal has_alcohol_tosylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                if reactants and product:
                    # Check for alcohol and tosyl chloride patterns in reactants
                    alcohol_pattern = Chem.MolFromSmarts("[OH][CH2]")
                    tosyl_chloride_pattern = Chem.MolFromSmarts("Cl[S](=[O])(=[O])[c]")

                    has_alcohol = False
                    has_tosyl_chloride = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(alcohol_pattern):
                                has_alcohol = True
                            if reactant_mol.HasSubstructMatch(tosyl_chloride_pattern):
                                has_tosyl_chloride = True

                    # Check for tosylate pattern in product
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        tosylate_pattern = Chem.MolFromSmarts("[O][S](=[O])(=[O])[c]")
                        if product_mol.HasSubstructMatch(tosylate_pattern):
                            if has_alcohol and has_tosyl_chloride:
                                has_alcohol_tosylation = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Alcohol tosylation strategy: {has_alcohol_tosylation}")
    return has_alcohol_tosylation
