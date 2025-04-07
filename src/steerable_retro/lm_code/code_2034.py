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
    Detects if the synthesis includes a Suzuki coupling between an aryl bromide and an aryl boronic ester.
    """
    found_suzuki = False

    def dfs_traverse(node, depth=0):
        nonlocal found_suzuki

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl bromide and boronic ester patterns in reactants
                aryl_bromide_found = False
                boronic_ester_found = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        aryl_bromide_pattern = Chem.MolFromSmarts("c[Br]")
                        boronic_ester_pattern = Chem.MolFromSmarts("cB(O)O")

                        if reactant_mol.HasSubstructMatch(aryl_bromide_pattern):
                            aryl_bromide_found = True
                        if reactant_mol.HasSubstructMatch(boronic_ester_pattern):
                            boronic_ester_found = True

                # If both patterns found, it's likely a Suzuki coupling
                if aryl_bromide_found and boronic_ester_found:
                    print(f"Found Suzuki coupling at depth {depth}")
                    found_suzuki = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_suzuki
