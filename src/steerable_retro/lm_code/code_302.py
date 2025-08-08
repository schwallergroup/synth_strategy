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
    Detects if the synthesis involves pyrazole formation in the early stages.
    Early stage is defined as high depth in the retrosynthetic tree.
    """
    pyrazole_formed = False
    high_depth_threshold = 3  # Consider depth >= 3 as early stage

    def dfs_traverse(node, depth=0):
        nonlocal pyrazole_formed

        if node["type"] == "reaction" and depth >= high_depth_threshold:
            # Check if this reaction forms a pyrazole
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains pyrazole but reactants don't
                product_mol = Chem.MolFromSmiles(product)
                pyrazole_pattern = Chem.MolFromSmarts("[n]1[n][c][c][c]1")

                if product_mol and product_mol.HasSubstructMatch(pyrazole_pattern):
                    # Check if reactants don't have pyrazole
                    has_pyrazole_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(pyrazole_pattern):
                            has_pyrazole_in_reactants = True
                            break

                    if not has_pyrazole_in_reactants:
                        print(f"Pyrazole formation detected at depth {depth}")
                        pyrazole_formed = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return pyrazole_formed
