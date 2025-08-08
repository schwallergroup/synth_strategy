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
    This function detects Suzuki coupling reactions in the synthetic route.
    Looks for aryl halide + boronic acid → biaryl transformations.
    """
    suzuki_found = False

    def dfs_traverse(node):
        nonlocal suzuki_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl halide (particularly bromide) in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")
            # Check for boronic acid in reactants
            boronic_acid_pattern = Chem.MolFromSmarts("[c]B(O)O")

            has_aryl_halide = False
            has_boronic_acid = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True
                    if mol and mol.HasSubstructMatch(boronic_acid_pattern):
                        has_boronic_acid = True
                except:
                    continue

            # Check if product has a biaryl bond that wasn't in the reactants
            if has_aryl_halide and has_boronic_acid:
                print("Potential Suzuki coupling detected")
                suzuki_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_found
