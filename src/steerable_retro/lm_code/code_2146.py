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
    Detects if the synthesis includes preparation of boronic acid for subsequent coupling.
    """
    boronic_acid_prep_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal boronic_acid_prep_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for aryl halide in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")

            # Check for boronic acid in product
            boronic_acid_pattern = Chem.MolFromSmarts("[c][B]([OH])[OH]")

            has_aryl_halide = False
            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                    has_aryl_halide = True
                    break

            product_mol = Chem.MolFromSmiles(product_smiles)
            has_boronic_acid = product_mol and product_mol.HasSubstructMatch(boronic_acid_pattern)

            if has_aryl_halide and has_boronic_acid:
                print(f"Boronic acid preparation detected at depth {depth}")
                boronic_acid_prep_detected = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return boronic_acid_prep_detected
