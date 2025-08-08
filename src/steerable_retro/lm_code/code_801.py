#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects if the synthesis uses a Suzuki coupling strategy
    to connect two heterocyclic fragments.
    """
    suzuki_found = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for boronic acid pattern
            boronic_acid_pattern = Chem.MolFromSmarts("[#6]B(O)O")

            # Check for aryl halide pattern
            aryl_halide_pattern = Chem.MolFromSmarts("[#6]~[#53,#35,#17]")

            # Check if reactants contain these patterns
            boronic_acid_present = False
            aryl_halide_present = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(boronic_acid_pattern):
                        boronic_acid_present = True
                    if mol.HasSubstructMatch(aryl_halide_pattern):
                        aryl_halide_present = True

            # If both patterns are present and product has a new C-C bond, likely a Suzuki coupling
            if boronic_acid_present and aryl_halide_present:
                print(f"Potential Suzuki coupling detected at depth {depth}")
                suzuki_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return suzuki_found
