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
    Detects heterocycle formation via diamine-dicarbonyl condensation,
    specifically looking for imidazopyridine formation.
    """
    found_pattern = False

    def dfs_traverse(node):
        nonlocal found_pattern

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for diamine pattern in reactants
            diamine_pattern = Chem.MolFromSmarts("[NH2][c][c][NH2]")
            dicarbonyl_pattern = Chem.MolFromSmarts("[CH]=[O].[CH]=[O]")

            # Check for imidazopyridine-like pattern in product
            imidazopyridine_pattern = Chem.MolFromSmarts("[#6]1[#7][#6][#6][#7][#6]1")

            reactant_has_diamine = False
            reactant_has_dicarbonyl = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(diamine_pattern):
                    reactant_has_diamine = True
                if mol and mol.HasSubstructMatch(dicarbonyl_pattern):
                    reactant_has_dicarbonyl = True

            product_mol = Chem.MolFromSmiles(product)
            product_has_imidazopyridine = product_mol and product_mol.HasSubstructMatch(
                imidazopyridine_pattern
            )

            if reactant_has_diamine and reactant_has_dicarbonyl and product_has_imidazopyridine:
                print("Found heterocycle formation via diamine-dicarbonyl condensation")
                found_pattern = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_pattern
