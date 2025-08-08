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
    This function detects a strategy involving aromatic halogenation followed by
    metal-catalyzed coupling (e.g., bromination followed by Suzuki coupling).
    """
    has_halogenation = False
    has_coupling = False

    def dfs_traverse(node):
        nonlocal has_halogenation, has_coupling

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for halogenation
            halide_pattern = Chem.MolFromSmarts("[c]-[Br,I,Cl]")

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(halide_pattern):
                # Check if reactants don't have the halide
                halide_in_reactants = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(halide_pattern):
                        halide_in_reactants = True
                        break

                if not halide_in_reactants:
                    has_halogenation = True
                    print(f"Detected aromatic halogenation: {rsmi}")

            # Check for coupling reaction
            boronic_pattern = Chem.MolFromSmarts("[c,C]-[B]([O])[O]")

            boronic_found = False
            halide_found = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(boronic_pattern):
                        boronic_found = True
                    if mol.HasSubstructMatch(halide_pattern):
                        halide_found = True

            if boronic_found and halide_found:
                has_coupling = True
                print(f"Detected coupling reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return has_halogenation and has_coupling
