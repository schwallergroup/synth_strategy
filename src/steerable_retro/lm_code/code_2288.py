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
    This function detects Suzuki coupling strategy where an aryl bromide
    reacts with a cyclopropylboronic acid.
    """
    suzuki_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactants contain aryl bromide and cyclopropylboronic acid
            aryl_bromide_pattern = Chem.MolFromSmarts("[c][Br]")
            boronic_acid_pattern = Chem.MolFromSmarts("[#6][B]")
            cyclopropyl_pattern = Chem.MolFromSmarts("[C]1[C][C]1")

            reactants_have_aryl_bromide = False
            reactants_have_cyclopropyl_boronic = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(aryl_bromide_pattern):
                        reactants_have_aryl_bromide = True
                    if (
                        mol
                        and mol.HasSubstructMatch(boronic_acid_pattern)
                        and mol.HasSubstructMatch(cyclopropyl_pattern)
                    ):
                        reactants_have_cyclopropyl_boronic = True
                except:
                    continue

            if reactants_have_aryl_bromide and reactants_have_cyclopropyl_boronic:
                print("Suzuki coupling with cyclopropylboronic acid detected")
                suzuki_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_detected
