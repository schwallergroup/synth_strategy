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
    This function detects cyano group introduction via aromatic substitution.
    """
    cyano_introduction_detected = False

    def dfs_traverse(node):
        nonlocal cyano_introduction_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl bromide to aryl cyano conversion
                aryl_bromide_pattern = Chem.MolFromSmarts("c-[Br]")
                aryl_cyano_pattern = Chem.MolFromSmarts("c-[C]#[N]")

                p_mol = Chem.MolFromSmiles(product)

                if p_mol and p_mol.HasSubstructMatch(aryl_cyano_pattern):
                    for reactant in reactants:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol and r_mol.HasSubstructMatch(aryl_bromide_pattern):
                            print(f"Cyano introduction detected: {rsmi}")
                            cyano_introduction_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return cyano_introduction_detected
