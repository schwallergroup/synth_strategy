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
    This function detects propargyl ether formation via alkylation of phenol.
    """
    propargyl_ether_found = False

    def dfs_traverse(node):
        nonlocal propargyl_ether_found

        if node["type"] == "reaction" and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phenol and propargyl bromide patterns in reactants
            phenol_pattern = Chem.MolFromSmarts("c[OH]")
            propargyl_bromide_pattern = Chem.MolFromSmarts("C#CC[Br]")
            propargyl_ether_pattern = Chem.MolFromSmarts("C#CCOc")

            phenol_found = False
            propargyl_bromide_found = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(phenol_pattern):
                    phenol_found = True
                if mol and mol.HasSubstructMatch(propargyl_bromide_pattern):
                    propargyl_bromide_found = True

            product_mol = Chem.MolFromSmiles(product)
            propargyl_ether_formed = product_mol and product_mol.HasSubstructMatch(
                propargyl_ether_pattern
            )

            if phenol_found and propargyl_bromide_found and propargyl_ether_formed:
                propargyl_ether_found = True
                print("Found propargyl ether formation via phenol alkylation")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return propargyl_ether_found
