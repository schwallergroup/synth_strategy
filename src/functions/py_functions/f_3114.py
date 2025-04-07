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
    This function detects if the synthesis route includes a Weinreb amide coupling
    to form an aryl ketone (aryl-C(=O)- bond formation).
    """
    weinreb_coupling_found = False

    def dfs_traverse(node):
        nonlocal weinreb_coupling_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl bromide pattern in reactants
            aryl_bromide_pattern = Chem.MolFromSmarts("[Br][c]")

            # Check for Weinreb amide pattern (simplified)
            weinreb_pattern = Chem.MolFromSmarts("[C](=O)N")

            # Check for aryl ketone pattern in product
            aryl_ketone_pattern = Chem.MolFromSmarts("[c][C](=O)[C]")

            reactant_has_aryl_bromide = False
            reactant_has_weinreb = False
            product_has_aryl_ketone = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(aryl_bromide_pattern):
                    reactant_has_aryl_bromide = True
                if mol and mol.HasSubstructMatch(weinreb_pattern):
                    reactant_has_weinreb = True

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(aryl_ketone_pattern):
                product_has_aryl_ketone = True

            if (
                reactant_has_aryl_bromide
                and reactant_has_weinreb
                and product_has_aryl_ketone
            ):
                print("Weinreb amide coupling detected")
                weinreb_coupling_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return weinreb_coupling_found
