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
    Detects if the synthetic route involves sequential functionalization
    of a thiophene core.
    """
    thiophene_core_found = False
    sequential_functionalization = False
    functionalization_steps = 0

    def dfs_traverse(node):
        nonlocal thiophene_core_found, sequential_functionalization, functionalization_steps

        if node["type"] == "mol":
            # Check if molecule contains thiophene core
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                thiophene_pattern = Chem.MolFromSmarts("[#6]1[#6][#16][#6][#6]1")
                if mol.HasSubstructMatch(thiophene_pattern):
                    thiophene_core_found = True

        elif node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactant_mol = Chem.MolFromSmiles(reactants_part)
            product_mol = Chem.MolFromSmiles(product_part)

            if reactant_mol and product_mol:
                thiophene_pattern = Chem.MolFromSmarts("[#6]1[#6][#16][#6][#6]1")

                # Check if both reactant and product contain thiophene
                if reactant_mol.HasSubstructMatch(
                    thiophene_pattern
                ) and product_mol.HasSubstructMatch(thiophene_pattern):
                    functionalization_steps += 1

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If we found thiophene core and at least 2 functionalization steps
    if thiophene_core_found and functionalization_steps >= 2:
        sequential_functionalization = True
        print(f"Found thiophene with {functionalization_steps} sequential functionalization steps")

    return sequential_functionalization
