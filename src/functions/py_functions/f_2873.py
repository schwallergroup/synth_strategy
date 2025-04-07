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
    Detects if the synthesis includes pyrazole ring formation via hydrazine-aldehyde condensation.
    """
    pyrazole_formation_detected = False

    def dfs_traverse(node):
        nonlocal pyrazole_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for hydrazine pattern in reactants
                hydrazine_pattern = Chem.MolFromSmarts("[#6]-[#7][#7]")
                aldehyde_pattern = Chem.MolFromSmarts("[#6][CH]=O")
                pyrazole_pattern = Chem.MolFromSmarts("[#6]1[n][n][c][c]1")

                has_hydrazine = any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(hydrazine_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                )
                has_aldehyde = any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(aldehyde_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                )

                product_mol = Chem.MolFromSmiles(product)
                has_pyrazole_in_product = product_mol and product_mol.HasSubstructMatch(
                    pyrazole_pattern
                )

                if has_hydrazine and has_aldehyde and has_pyrazole_in_product:
                    print(
                        "Pyrazole formation via hydrazine-aldehyde condensation detected"
                    )
                    pyrazole_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Pyrazole formation via condensation: {pyrazole_formation_detected}")
    return pyrazole_formation_detected
