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
    Detects synthesis where fragments containing heterocycles (oxazole, pyrazole)
    are coupled in a late stage.
    """
    # Track if we found the pattern
    found_pattern = False
    # Track if fragments contain heterocycles
    fragment_has_oxazole = False
    fragment_has_pyrazole = False
    # Track if coupling occurs late stage
    late_stage_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern, fragment_has_oxazole, fragment_has_pyrazole, late_stage_coupling

        if node["type"] == "mol":
            # Check molecule for heterocycles
            if "smiles" in node:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for oxazole
                    oxazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#6][#6][#8]1")
                    if mol.HasSubstructMatch(oxazole_pattern):
                        fragment_has_oxazole = True
                        print(f"Found oxazole-containing fragment")

                    # Check for pyrazole
                    pyrazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#7][#6][#6]1")
                    if mol.HasSubstructMatch(pyrazole_pattern):
                        fragment_has_pyrazole = True
                        print(f"Found pyrazole-containing fragment")

        elif node["type"] == "reaction" and depth <= 1:  # Late stage reaction
            # Check if this is a coupling reaction
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                if len(reactants) > 1:
                    late_stage_coupling = True
                    print(f"Found late-stage coupling at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have the pattern we're looking for
    if fragment_has_oxazole and fragment_has_pyrazole and late_stage_coupling:
        found_pattern = True
        print("Found coupling of heterocycle-containing fragments")

    return found_pattern
