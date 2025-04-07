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
    Detects if the synthesis includes strategic alcohol functional group transformations
    (alcohol → aldehyde and/or alcohol → bromide)
    """
    alcohol_to_aldehyde = False
    alcohol_to_bromide = False

    def dfs_traverse(node):
        nonlocal alcohol_to_aldehyde, alcohol_to_bromide

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                product_mol = Chem.MolFromSmiles(product_part)

                if product_mol:
                    # Patterns
                    alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8H]")
                    aldehyde_pattern = Chem.MolFromSmarts("[#6H]-[#8]=[#6]")
                    bromide_pattern = Chem.MolFromSmarts("[#6]-[#35]")

                    # Check for alcohol to aldehyde
                    if any(
                        r and r.HasSubstructMatch(alcohol_pattern)
                        for r in reactants
                        if r
                    ):
                        if product_mol.HasSubstructMatch(aldehyde_pattern):
                            print("Detected alcohol to aldehyde transformation")
                            alcohol_to_aldehyde = True

                    # Check for alcohol to bromide
                    if any(
                        r and r.HasSubstructMatch(alcohol_pattern)
                        for r in reactants
                        if r
                    ):
                        if product_mol.HasSubstructMatch(bromide_pattern):
                            print("Detected alcohol to bromide transformation")
                            alcohol_to_bromide = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    # Return True if either transformation is detected
    return alcohol_to_aldehyde or alcohol_to_bromide
