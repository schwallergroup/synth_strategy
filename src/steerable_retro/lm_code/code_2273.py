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
    This function detects a synthetic sequence involving nitrile and oxime transformations.
    """
    nitrile_formed = False
    nitrile_consumed = False
    oxime_formed = False
    oxime_consumed = False

    def dfs_traverse(node):
        nonlocal nitrile_formed, nitrile_consumed, oxime_formed, oxime_consumed

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactant_mol = Chem.MolFromSmiles(reactants)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    nitrile_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
                    oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2]-[OX2]")

                    # Check for nitrile formation
                    if not reactant_mol.HasSubstructMatch(
                        nitrile_pattern
                    ) and product_mol.HasSubstructMatch(nitrile_pattern):
                        nitrile_formed = True
                        print("Detected nitrile formation")

                    # Check for nitrile consumption
                    if reactant_mol.HasSubstructMatch(
                        nitrile_pattern
                    ) and not product_mol.HasSubstructMatch(nitrile_pattern):
                        nitrile_consumed = True
                        print("Detected nitrile consumption")

                    # Check for oxime formation
                    if not reactant_mol.HasSubstructMatch(
                        oxime_pattern
                    ) and product_mol.HasSubstructMatch(oxime_pattern):
                        oxime_formed = True
                        print("Detected oxime formation")

                    # Check for oxime consumption
                    if reactant_mol.HasSubstructMatch(
                        oxime_pattern
                    ) and not product_mol.HasSubstructMatch(oxime_pattern):
                        oxime_consumed = True
                        print("Detected oxime consumption")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we have both nitrile and oxime transformations
    return nitrile_formed and nitrile_consumed and oxime_formed and oxime_consumed
