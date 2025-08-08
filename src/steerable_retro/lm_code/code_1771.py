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
    Detects synthesis strategy involving activation of benzyl alcohol to benzyl halide
    for subsequent coupling reactions.
    """
    benzyl_alcohol_to_halide = False

    def dfs_traverse(node):
        nonlocal benzyl_alcohol_to_halide

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for benzyl alcohol to benzyl halide transformation
            benzyl_alcohol_pattern = Chem.MolFromSmarts("[c:1]-[CH2:2]-[OH:3]")
            benzyl_halide_pattern = Chem.MolFromSmarts("[c:1]-[CH2:2]-[Br,Cl,I:3]")

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(benzyl_alcohol_pattern):
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(benzyl_halide_pattern):
                        print("Detected benzyl alcohol to halide activation")
                        benzyl_alcohol_to_halide = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return benzyl_alcohol_to_halide
