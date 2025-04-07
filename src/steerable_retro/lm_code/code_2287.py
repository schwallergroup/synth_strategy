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
    This function detects SNAr diaryl ether formation strategy where an aryl fluoride
    reacts with a phenol to form a diaryl ether.
    """
    snar_detected = False

    def dfs_traverse(node):
        nonlocal snar_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactants contain aryl fluoride and phenol
            aryl_fluoride_pattern = Chem.MolFromSmarts("[c][F]")
            phenol_pattern = Chem.MolFromSmarts("[c][OH]")
            diaryl_ether_pattern = Chem.MolFromSmarts("[c][O][c]")

            reactants_have_aryl_fluoride = False
            reactants_have_phenol = False
            product_has_diaryl_ether = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(aryl_fluoride_pattern):
                        reactants_have_aryl_fluoride = True
                    if mol and mol.HasSubstructMatch(phenol_pattern):
                        reactants_have_phenol = True
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(diaryl_ether_pattern):
                    product_has_diaryl_ether = True
            except:
                pass

            if reactants_have_aryl_fluoride and reactants_have_phenol and product_has_diaryl_ether:
                print("SNAr diaryl ether formation detected")
                snar_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return snar_detected
