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
    This function detects Suzuki coupling reactions in the synthetic route.
    Looks for aryl-aryl bond formation between aryl halide and boronic acid.
    """
    suzuki_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid pattern in reactants
                boronic_acid_pattern = Chem.MolFromSmarts("[c;H0]-[B]([O])([O])")
                aryl_halide_pattern = Chem.MolFromSmarts("[c;H0]-[Cl,Br,I]")

                has_boronic_acid = False
                has_aryl_halide = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(boronic_acid_pattern):
                            has_boronic_acid = True
                        if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True
                    except:
                        continue

                # If both patterns are found in reactants, check for biaryl in product
                if has_boronic_acid and has_aryl_halide:
                    try:
                        prod_mol = Chem.MolFromSmiles(product)
                        # Biaryl formation is likely if Suzuki conditions are present
                        if prod_mol:
                            print("Potential Suzuki coupling detected")
                            suzuki_detected = True
                    except:
                        pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_detected
