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
    This function detects the esterification of a carboxylic acid in the synthesis route.
    """
    esterification_found = False

    def dfs_traverse(node):
        nonlocal esterification_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carboxylic acid in reactants
                carboxylic_acid_pattern = Chem.MolFromSmarts("[#6]-C(=O)[OH]")
                alcohol_pattern = Chem.MolFromSmarts("[#6]-[OH]")
                ester_pattern = Chem.MolFromSmarts("[#6]-[#8]-C(=O)-[#6]")

                has_acid = False
                has_alcohol = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(carboxylic_acid_pattern):
                                has_acid = True
                            if mol.HasSubstructMatch(alcohol_pattern):
                                has_alcohol = True
                    except:
                        continue

                # Check if product has ester
                if has_acid and has_alcohol:
                    try:
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol and prod_mol.HasSubstructMatch(ester_pattern):
                            esterification_found = True
                            print("Found esterification reaction")
                    except:
                        pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return esterification_found
