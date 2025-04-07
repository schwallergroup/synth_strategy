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
    This function detects a sequence involving alcohol protection with TBDPS
    and oxidation/reduction steps.
    """
    protection_step = False
    oxidation_step = False

    def dfs_traverse(node):
        nonlocal protection_step, oxidation_step

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for TBDPS protection
                product_mol = Chem.MolFromSmiles(product)
                tbdps_pattern = Chem.MolFromSmarts(
                    "[O][Si]([c]1[cH][cH][cH][cH][cH]1)([c]1[cH][cH][cH][cH][cH]1)[C]([CH3])([CH3])[CH3]"
                )
                alcohol_pattern = Chem.MolFromSmarts("[OH][CH]")

                if product_mol and product_mol.HasSubstructMatch(tbdps_pattern):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            alcohol_pattern
                        ):
                            protection_step = True
                            print("TBDPS protection detected")

                # Check for alcohol oxidation
                ketone_pattern = Chem.MolFromSmarts("[C](=[O])[C]")
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(ketone_pattern):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            alcohol_pattern
                        ):
                            oxidation_step = True
                            print("Alcohol oxidation detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return protection_step and oxidation_step
