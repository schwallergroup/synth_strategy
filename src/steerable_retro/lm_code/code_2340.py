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
    This function detects if the route contains aromatization of a saturated heterocycle,
    specifically tetrahydropyridine to pyridine conversion.
    """
    found = False

    def dfs_traverse(node):
        nonlocal found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for tetrahydropyridine pattern in reactants
            reactant_mol = Chem.MolFromSmiles(reactants)
            tetrahydropyridine_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#7][#6][#6]1")

            # Check for pyridine pattern in product
            product_mol = Chem.MolFromSmiles(product)
            pyridine_pattern = Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#7]:[#6]:[#6]:1")

            if reactant_mol and product_mol:
                if reactant_mol.HasSubstructMatch(
                    tetrahydropyridine_pattern
                ) and product_mol.HasSubstructMatch(pyridine_pattern):
                    print("Found heterocycle aromatization (tetrahydropyridine to pyridine)")
                    found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return found
