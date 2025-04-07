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
    This function detects formation of sulfonyl chloride from sulfonic acid.
    """
    sulfonyl_chloride_formation_detected = False

    def dfs_traverse(node):
        nonlocal sulfonyl_chloride_formation_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Try to detect sulfonyl chloride formation pattern
            try:
                p_mol = Chem.MolFromSmiles(product)
                sulfonyl_chloride_pattern = Chem.MolFromSmarts("[S](=O)(=O)[Cl]")

                if p_mol and p_mol.HasSubstructMatch(sulfonyl_chloride_pattern):
                    for r in reactants:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol:
                            sulfonic_acid_pattern = Chem.MolFromSmarts("[S](=O)(=O)[O]")
                            if r_mol.HasSubstructMatch(sulfonic_acid_pattern):
                                print(
                                    "Sulfonyl chloride formation detected in reaction:",
                                    rsmi,
                                )
                                sulfonyl_chloride_formation_detected = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return sulfonyl_chloride_formation_detected
