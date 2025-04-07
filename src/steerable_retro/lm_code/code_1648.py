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
    This function detects the use of nitrile-containing heterocycles in the synthesis.
    """
    has_nitrile_heterocycle = False

    def dfs_traverse(node):
        nonlocal has_nitrile_heterocycle

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check for nitrile-containing heterocycle in reactants
                nitrile_pattern = Chem.MolFromSmarts("C#N")
                heterocycle_patterns = [
                    Chem.MolFromSmarts("c1nccn1"),  # pyrazole
                    Chem.MolFromSmarts("c1ncnn1"),  # triazole
                    Chem.MolFromSmarts("c1ncccc1"),  # pyridine
                    Chem.MolFromSmarts("c1nccs1"),  # thiazole
                    Chem.MolFromSmarts("c1nccnc1"),  # pyrimidine
                ]

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(nitrile_pattern):
                            # Check if it also has a heterocycle
                            for pattern in heterocycle_patterns:
                                if mol.HasSubstructMatch(pattern):
                                    print("Detected nitrile-containing heterocycle")
                                    has_nitrile_heterocycle = True
                                    break
                            if has_nitrile_heterocycle:
                                break
                    except:
                        continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return has_nitrile_heterocycle
