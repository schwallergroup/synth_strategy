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
    This function detects diaryl ether formation via SNAr reaction.
    Looks for reactions where a halogenated aromatic compound reacts with a phenol.
    """
    snar_found = False

    def dfs_traverse(node):
        nonlocal snar_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")

            # Check if one reactant has aromatic halide and another has phenol
            has_aryl_halide = False
            has_phenol = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    # Check for aromatic halide
                    aryl_halide_pattern = Chem.MolFromSmarts("[c][Cl,F,Br,I]")
                    if mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True

                    # Check for phenol
                    phenol_pattern = Chem.MolFromSmarts("[c][OH]")
                    if mol.HasSubstructMatch(phenol_pattern):
                        has_phenol = True

            # Check if product has diaryl ether
            product_mol = Chem.MolFromSmiles(product_part)
            if product_mol:
                diaryl_ether_pattern = Chem.MolFromSmarts("[c][O][c]")
                if product_mol.HasSubstructMatch(diaryl_ether_pattern):
                    # If reactants had aryl halide and phenol, and product has diaryl ether, it's likely SNAr
                    if has_aryl_halide and has_phenol:
                        print("Found diaryl ether formation via SNAr")
                        snar_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return snar_found
