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
    Detects if the synthesis follows a linear strategy (each step adds one fragment)
    while maintaining a pyridine core throughout.
    """
    maintains_pyridine_core = True
    is_linear_synthesis = True

    def dfs_traverse(node):
        nonlocal maintains_pyridine_core, is_linear_synthesis

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product has pyridine core
            product_mol = Chem.MolFromSmiles(product)
            pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")

            if product_mol and not product_mol.HasSubstructMatch(pyridine_pattern):
                maintains_pyridine_core = False

            # Check if synthesis is linear (no more than 2 significant fragments)
            significant_fragments = 0
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if (
                    mol and mol.GetNumHeavyAtoms() > 3
                ):  # Consider fragments with >3 atoms as significant
                    significant_fragments += 1

            if significant_fragments > 2:
                is_linear_synthesis = False

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return maintains_pyridine_core and is_linear_synthesis
