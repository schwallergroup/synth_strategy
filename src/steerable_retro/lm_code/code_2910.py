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
    Detects if the synthesis route uses halogen displacement chemistry
    """
    halogen_displacement_count = 0

    def dfs_traverse(node):
        nonlocal halogen_displacement_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactants contain halogens
            halogen_in_reactants = False
            for reactant in reactants:
                if "Br" in reactant or "Cl" in reactant or "I" in reactant:
                    halogen_in_reactants = True
                    break

            # Simplified check - in a real implementation, you would need to verify
            # that the halogen is actually being displaced
            if halogen_in_reactants and ("[N:" in rsmi or "[C:" in rsmi):
                print("Potential halogen displacement detected")
                halogen_displacement_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return halogen_displacement_count >= 1
