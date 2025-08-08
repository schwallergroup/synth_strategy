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
    This function detects a synthetic strategy involving reaction between a thiol and isocyanate.
    """
    thiol_isocyanate_detected = False

    def dfs_traverse(node):
        nonlocal thiol_isocyanate_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if reactants contain thiol and isocyanate
                thiol_pattern = Chem.MolFromSmarts("[#16H]")
                isocyanate_pattern = Chem.MolFromSmarts("[#7]=[#6]=[#8]")

                has_thiol = False
                has_isocyanate = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol is not None:
                        if mol.HasSubstructMatch(thiol_pattern):
                            has_thiol = True
                        if mol.HasSubstructMatch(isocyanate_pattern):
                            has_isocyanate = True

                if has_thiol and has_isocyanate:
                    print("Found thiol-isocyanate reaction")
                    thiol_isocyanate_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thiol_isocyanate_detected
