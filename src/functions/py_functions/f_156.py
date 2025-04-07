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
    Detects if the synthesis includes formation of a morpholine amide
    """
    found_morpholine_amide_formation = False

    def dfs_traverse(node):
        nonlocal found_morpholine_amide_formation

        if (
            node.get("type") == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for morpholine amide formation
            morpholine_pattern = Chem.MolFromSmarts("N1CCOCC1")
            morpholine_amide_pattern = Chem.MolFromSmarts("C(=O)N1CCOCC1")
            acid_chloride_pattern = Chem.MolFromSmarts("C(=O)Cl")

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(morpholine_amide_pattern):
                has_morpholine = False
                has_acid_chloride = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(morpholine_pattern):
                            has_morpholine = True
                        if reactant_mol.HasSubstructMatch(acid_chloride_pattern):
                            has_acid_chloride = True

                if has_morpholine and has_acid_chloride:
                    found_morpholine_amide_formation = True
                    print("Detected morpholine amide formation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_morpholine_amide_formation
