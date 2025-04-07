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
    This function detects if the synthesis includes a C-C bond formation between two heterocyclic systems.
    """
    heterocycle_coupling_detected = False

    def dfs_traverse(node):
        nonlocal heterocycle_coupling_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if we have at least two reactants (potential for coupling)
            if len(reactants) >= 2:
                # Check if both reactants contain heterocycles
                heterocycle_patterns = [
                    Chem.MolFromSmarts("[#6]1:[#7]:[#6]:[#6]:[#6]:[#6]:1"),  # pyridine
                    Chem.MolFromSmarts("[#6]1:[#7]:[#6]:[#7]:[#6]:1"),  # pyrimidine
                    Chem.MolFromSmarts("[#6]1:[#7]:[#6]:[#6]:[#7]:1"),  # pyrazine
                    Chem.MolFromSmarts("[#6]1:[#7]:[#7]:[#6]:[#6]:1"),  # pyridazine
                    Chem.MolFromSmarts(
                        "[#6]1:[#7]:[#6]:[#6]:[#6]:[#7]:1"
                    ),  # quinoxaline
                ]

                heterocycle_count = 0
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    for pattern in heterocycle_patterns:
                        if reactant_mol.HasSubstructMatch(pattern):
                            heterocycle_count += 1
                            break

                product_mol = Chem.MolFromSmiles(product)
                if heterocycle_count >= 2 and product_mol:
                    # Check if product has a more complex heterocyclic system
                    for pattern in heterocycle_patterns:
                        if product_mol.HasSubstructMatch(pattern):
                            heterocycle_coupling_detected = True
                            print("Detected heterocycle coupling")
                            break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return heterocycle_coupling_detected
