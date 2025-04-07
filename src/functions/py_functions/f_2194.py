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
    Detects if the synthesis includes formation of a pyrazole ring
    from a dicarbonyl compound and hydrazine.
    """
    found_pyrazole_formation = False

    def dfs_traverse(node):
        nonlocal found_pyrazole_formation

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for hydrazine pattern in reactants
            hydrazine_pattern = Chem.MolFromSmarts("[N][N]")
            dicarbonyl_pattern = Chem.MolFromSmarts("[C](=[O])[C][C](=[O])")

            # Check for pyrazole pattern in product
            pyrazole_pattern = Chem.MolFromSmarts("[n]1[n][c][c][c]1")

            has_hydrazine = False
            has_dicarbonyl = False
            forms_pyrazole = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(hydrazine_pattern):
                        has_hydrazine = True
                    if mol and mol.HasSubstructMatch(dicarbonyl_pattern):
                        has_dicarbonyl = True
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(pyrazole_pattern):
                    forms_pyrazole = True
            except:
                pass

            if has_hydrazine and has_dicarbonyl and forms_pyrazole:
                print("Found pyrazole formation from hydrazine and dicarbonyl")
                found_pyrazole_formation = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if found_pyrazole_formation:
        print("Detected pyrazole formation strategy")

    return found_pyrazole_formation
