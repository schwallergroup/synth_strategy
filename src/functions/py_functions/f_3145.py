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
    Detects pyrazole ring formation using hydrazine as a reagent.
    Looks for a reaction where hydrazine (NH2-NH2) is used to form a 5-membered heterocyclic ring.
    """
    pyrazole_formation_detected = False

    def dfs_traverse(node):
        nonlocal pyrazole_formation_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if hydrazine is a reactant
                hydrazine_pattern = Chem.MolFromSmarts("[NH2][NH2]")
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(hydrazine_pattern):
                            # Check if product has a pyrazole ring
                            product_mol = Chem.MolFromSmiles(product)
                            pyrazole_pattern = Chem.MolFromSmarts("c1nn[c,n]c1")
                            if product_mol and product_mol.HasSubstructMatch(
                                pyrazole_pattern
                            ):
                                print("Detected pyrazole formation using hydrazine")
                                pyrazole_formation_detected = True
                    except:
                        continue

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return pyrazole_formation_detected
