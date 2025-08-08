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
    Detects if the synthesis includes N-alkylation with a benzyl halide.
    Looks for [NH]-heterocycle + benzyl-X → N-benzyl-heterocycle pattern.
    """
    found = False

    def dfs_traverse(node):
        nonlocal found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for NH-heterocycle pattern
            nh_pattern = Chem.MolFromSmarts("[nH]")
            # Check for benzyl halide pattern
            benzyl_halide_pattern = Chem.MolFromSmarts("c[CH2][Cl,Br,I]")
            # Check for N-benzyl pattern in product
            n_benzyl_pattern = Chem.MolFromSmarts("n[CH2]c")

            has_nh = False
            has_benzyl_halide = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(nh_pattern):
                        has_nh = True
                    if mol.HasSubstructMatch(benzyl_halide_pattern):
                        has_benzyl_halide = True

            if has_nh and has_benzyl_halide:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and prod_mol.HasSubstructMatch(n_benzyl_pattern):
                    print("Found N-alkylation with benzyl halide")
                    found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found
