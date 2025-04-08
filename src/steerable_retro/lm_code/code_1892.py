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
    Detects if the synthetic route involves tetrazole ring formation.
    """
    tetrazole_formation_found = False

    def dfs_traverse(node):
        nonlocal tetrazole_formation_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains tetrazole but reactants don't
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol:
                tetrazole_pattern = Chem.MolFromSmarts("[#7]1[#7][#7][#7][#6]1")

                if product_mol.HasSubstructMatch(tetrazole_pattern):
                    # Check if reactants don't have tetrazole
                    tetrazole_in_reactants = False
                    for r in reactants:
                        reactant_mol = Chem.MolFromSmiles(r) if r else None
                        if reactant_mol and reactant_mol.HasSubstructMatch(tetrazole_pattern):
                            tetrazole_in_reactants = True
                            break

                    if not tetrazole_in_reactants:
                        print("Found tetrazole formation")
                        tetrazole_formation_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return tetrazole_formation_found
