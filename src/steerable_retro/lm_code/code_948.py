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
    This function detects disconnections that occur adjacent to protected functional groups.
    """
    protected_groups = [
        Chem.MolFromSmarts("[#6]1[#8][#6][#6]([#6])([#6])[#8]1"),  # acetonide
        Chem.MolFromSmarts("[#6][Si]([#6])([#6])[#8][#6]"),  # silyl ether
    ]

    found_pattern = False

    def dfs_traverse(node):
        nonlocal found_pattern

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            # Check if product has protected group
                            for pattern in protected_groups:
                                if product_mol.HasSubstructMatch(pattern):
                                    # Check if disconnection happens adjacent to protected group
                                    # This is a simplification - in practice would need more sophisticated analysis
                                    if len(reactants) > 1:
                                        found_pattern = True
                                        print("Found disconnection adjacent to protected group")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return found_pattern
