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
    Detects if the synthesis route involves borylation of an aryl halide
    to prepare for a subsequent coupling reaction.
    """
    borylation_detected = False

    def dfs_traverse(node):
        nonlocal borylation_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for borylation pattern
            product_mol = Chem.MolFromSmiles(product)
            boronate_patt = Chem.MolFromSmarts("[#6]-[B]([O][C])[O][C]")

            if product_mol and product_mol.HasSubstructMatch(boronate_patt):
                # Check if reactant had a halide
                halide_patt = Chem.MolFromSmarts("[#6]-[Br,I,Cl]")
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(halide_patt):
                        print("Borylation detected")
                        borylation_detected = True
                        break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return borylation_detected
