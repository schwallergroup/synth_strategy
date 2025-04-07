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
    This function detects pyrazole ring formation from nitrile and hydrazine components.
    """
    pyrazole_formation_detected = False

    def dfs_traverse(node):
        nonlocal pyrazole_formation_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains pyrazole
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    pyrazole_pattern = Chem.MolFromSmarts("[n]1[n][c][c][c]1")
                    if product_mol.HasSubstructMatch(pyrazole_pattern):

                        # Check if reactants contain nitrile and hydrazine/hydrazide
                        nitrile_found = False
                        hydrazine_found = False

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                                hydrazine_pattern = Chem.MolFromSmarts("[NH]-[NH2]")

                                if reactant_mol.HasSubstructMatch(nitrile_pattern):
                                    nitrile_found = True
                                if reactant_mol.HasSubstructMatch(hydrazine_pattern):
                                    hydrazine_found = True

                        if nitrile_found and hydrazine_found:
                            pyrazole_formation_detected = True
                            print("Detected pyrazole formation via nitrile cyclization")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return pyrazole_formation_detected
