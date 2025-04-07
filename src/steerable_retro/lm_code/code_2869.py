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
    Detects if the synthetic route includes amide formation using an acid chloride.
    Looks for reaction between acid chloride and amine to form amide.
    """
    found_amide_formation = False

    def is_amide_formation(reactants, product):
        # Check if one reactant has acid chloride group
        acid_chloride_pattern = Chem.MolFromSmarts("[#6]-[C](=[O])[Cl]")
        has_acid_chloride = False
        for r in reactants:
            r_mol = Chem.MolFromSmiles(r)
            if r_mol and r_mol.HasSubstructMatch(acid_chloride_pattern):
                has_acid_chloride = True
                break

        if not has_acid_chloride:
            return False

        # Check if another reactant has amine group
        amine_pattern = Chem.MolFromSmarts("[NH2]")
        has_amine = False
        for r in reactants:
            r_mol = Chem.MolFromSmiles(r)
            if r_mol and r_mol.HasSubstructMatch(amine_pattern):
                has_amine = True
                break

        if not has_amine:
            return False

        # Check if product has amide group
        amide_pattern = Chem.MolFromSmarts("[#6]-[C](=[O])[NH]-[#6]")
        product_mol = Chem.MolFromSmiles(product)
        if product_mol and product_mol.HasSubstructMatch(amide_pattern):
            return True

        return False

    def dfs_traverse(node):
        nonlocal found_amide_formation

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            if is_amide_formation(reactants, product):
                found_amide_formation = True
                print("Found amide formation from acid chloride")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Amide formation from acid chloride: {found_amide_formation}")
    return found_amide_formation
