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
    This function detects Suzuki coupling reactions forming biaryl systems.
    Looks for reactions where a boronic acid and aryl halide form a biaryl system.
    """
    suzuki_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for boronic acid in reactants
            boronic_acid_pattern = Chem.MolFromSmarts("[c][B]([OH])[OH]")
            # Check for aryl halide in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")
            # Check for biaryl in product
            biaryl_pattern = Chem.MolFromSmarts("[c]!@[c]")

            has_boronic_acid = False
            has_aryl_halide = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(boronic_acid_pattern):
                        has_boronic_acid = True
                    if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                has_biaryl = product_mol and product_mol.HasSubstructMatch(biaryl_pattern)
            except:
                has_biaryl = False

            if has_boronic_acid and has_aryl_halide and has_biaryl:
                print("Suzuki coupling detected: boronic acid + aryl halide → biaryl")
                suzuki_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_detected
