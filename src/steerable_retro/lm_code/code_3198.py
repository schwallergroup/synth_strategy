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
    This function detects if the synthesis involves a carbonyl reduction (C=O to CH2-OH)
    in the early stages of the synthesis.
    """
    found_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_reduction

        if node["type"] == "reaction" and depth >= 2:  # Early stage (higher depth)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains C=O and any reactant contains CH2-OH
                product_mol = Chem.MolFromSmiles(product)
                carbonyl_pattern = Chem.MolFromSmarts("[#6]=O")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    alcohol_pattern = Chem.MolFromSmarts("[#6][OH]")

                    if (
                        product_mol
                        and reactant_mol
                        and product_mol.HasSubstructMatch(carbonyl_pattern)
                        and reactant_mol.HasSubstructMatch(alcohol_pattern)
                    ):
                        print(f"Found carbonyl reduction at depth {depth}")
                        found_reduction = True
                        break

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_reduction
