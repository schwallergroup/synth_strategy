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
    Detects if the synthesis uses a cross-coupling reaction (like Suzuki coupling).
    Looks for reactions where a boronic acid and aryl halide form a biaryl system.
    """
    found_cross_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_cross_coupling

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Patterns for boronic acid and aryl halide
                boronic_acid_pattern = Chem.MolFromSmarts("[#6][B]([OH])[OH]")
                aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")

                # Check for boronic acid and aryl halide in reactants
                has_boronic_acid = False
                has_aryl_halide = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if not mol:
                        continue
                    if mol.HasSubstructMatch(boronic_acid_pattern):
                        has_boronic_acid = True
                    if mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True

                # Check for biaryl in product
                product_mol = Chem.MolFromSmiles(product)
                biaryl_pattern = Chem.MolFromSmarts("c:c-c:c")  # Simple biaryl pattern
                has_biaryl = product_mol and product_mol.HasSubstructMatch(biaryl_pattern)

                if has_boronic_acid and has_aryl_halide and has_biaryl:
                    found_cross_coupling = True
                    print(f"Detected cross-coupling reaction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_cross_coupling
