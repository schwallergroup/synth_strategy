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
    This function detects if the synthetic route involves benzyl ether protection
    of a phenol in the early stages of synthesis.
    """
    benzyl_protection_found = False
    phenol_present = False

    def dfs_traverse(node, depth=0):
        nonlocal benzyl_protection_found, phenol_present

        if node["type"] == "reaction" and depth >= 2:  # Early in synthesis (higher depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this reaction involves benzyl protection
                benzyl_pattern = Chem.MolFromSmarts("[c][O][C][c]1[cH][cH][cH][cH][cH]1")
                phenol_pattern = Chem.MolFromSmarts("[c][OH]")

                product_mol = Chem.MolFromSmiles(product)

                # Check if any reactant has a phenol
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(phenol_pattern):
                        phenol_present = True

                # Check if product has benzyl ether
                if product_mol and product_mol.HasSubstructMatch(benzyl_pattern):
                    # If a reactant had phenol and product has benzyl ether, it's a protection
                    if phenol_present:
                        benzyl_protection_found = True
                        print(f"Benzyl protection found at depth {depth}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return benzyl_protection_found
