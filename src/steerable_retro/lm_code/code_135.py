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
    This function detects a synthetic strategy involving coupling with an isocyanate group.
    """
    isocyanate_coupling_detected = False

    def dfs_traverse(node):
        nonlocal isocyanate_coupling_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if any reactant contains isocyanate group
                isocyanate_pattern = Chem.MolFromSmarts("[#7]=[#6]=[#8]")

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol is not None and mol.HasSubstructMatch(isocyanate_pattern):
                        print("Found isocyanate coupling")
                        isocyanate_coupling_detected = True
                        break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return isocyanate_coupling_detected
