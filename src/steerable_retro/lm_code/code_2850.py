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
    Detects Suzuki coupling between heteroaromatic systems.
    Looks for C-C bond formation between a boronate and an aryl halide,
    where at least one is a heterocycle.
    """
    suzuki_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronate in reactants
                boronate_pattern = Chem.MolFromSmarts("[c][B][O][C]")
                # Check for aryl halide in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")
                # Check for heterocycles
                thiophene_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#16][#6]1")
                pyridine_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#7][#6]1")

                has_boronate = False
                has_aryl_halide = False
                has_heterocycle = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(boronate_pattern):
                            has_boronate = True
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True
                        if mol.HasSubstructMatch(thiophene_pattern) or mol.HasSubstructMatch(
                            pyridine_pattern
                        ):
                            has_heterocycle = True

                # Check if product has a new C-C bond between aromatics
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and has_boronate and has_aryl_halide and has_heterocycle:
                    # This is a simplified check - in a real implementation,
                    # we would need to track the specific atoms involved in the new bond
                    suzuki_detected = True
                    print(f"Detected Suzuki coupling between heteroaromatics: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_detected
