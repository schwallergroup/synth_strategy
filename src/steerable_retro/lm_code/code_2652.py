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
    This function detects if the route contains at least two Suzuki coupling reactions.
    Suzuki couplings typically involve an aryl halide and a boronic acid/ester.
    """
    suzuki_count = 0

    def dfs_traverse(node):
        nonlocal suzuki_count

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for boronic acid/ester pattern in reactants
            has_boron = False
            has_halide = False

            for reactant in reactants:
                if reactant:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Check for boronic acid/ester
                        if mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")
                        ) or mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]-[#5](-[#8;R])-[#8;R]")):
                            has_boron = True

                        # Check for aryl halide
                        if mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[#6;a]-[#35,#53]")
                        ) or mol.HasSubstructMatch(Chem.MolFromSmarts("[#6;a]-[#17]")):
                            has_halide = True

            # If both patterns are found, it's likely a Suzuki coupling
            if has_boron and has_halide:
                suzuki_count += 1
                print(f"Found Suzuki coupling: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Total Suzuki couplings found: {suzuki_count}")
    return suzuki_count >= 2
