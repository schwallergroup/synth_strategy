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
    Detects if the synthesis involves formation of an acid chloride from a carboxylic acid.
    """
    # Track if we found acid chloride formation
    found_acid_chloride_formation = False

    def dfs_traverse(node):
        nonlocal found_acid_chloride_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if product contains acid chloride
            product_mol = Chem.MolFromSmiles(product_part)
            if product_mol:
                acid_chloride_pattern = Chem.MolFromSmarts("[#6](=[#8])[Cl]")
                if product_mol.HasSubstructMatch(acid_chloride_pattern):
                    # Check if reactants contain carboxylic acid
                    reactants = reactants_part.split(".")
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            carboxylic_acid_pattern = Chem.MolFromSmarts("[#6](=[#8])[#8H]")
                            if reactant_mol.HasSubstructMatch(carboxylic_acid_pattern):
                                print("Found acid chloride formation from carboxylic acid")
                                found_acid_chloride_formation = True
                                break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return found_acid_chloride_formation
