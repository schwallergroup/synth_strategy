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
    This function detects if the synthetic route employs a Suzuki coupling for biaryl formation.
    """
    suzuki_found = False

    def dfs_traverse(node):
        nonlocal suzuki_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for boronic acid in reactants
            boronic_acid_pattern = Chem.MolFromSmarts("[c][B]([OH])[OH]")
            # Check for aryl halide in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("[c][Cl,Br,I]")
            # Check for biaryl in product
            biaryl_pattern = Chem.MolFromSmarts("[c](-[c])[c]")

            boronic_acid_present = False
            aryl_halide_present = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(boronic_acid_pattern):
                        boronic_acid_present = True
                    if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                        aryl_halide_present = True
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                biaryl_present = product_mol and product_mol.HasSubstructMatch(biaryl_pattern)
            except:
                biaryl_present = False

            if boronic_acid_present and aryl_halide_present and biaryl_present:
                print("Suzuki coupling detected")
                suzuki_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_found
