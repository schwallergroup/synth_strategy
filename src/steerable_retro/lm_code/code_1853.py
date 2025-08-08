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
    Detects a strategy involving a Suzuki coupling in the second half of the synthesis.
    """
    # Track if we find Suzuki coupling reactants and products
    found_suzuki = False
    suzuki_depth = None
    max_depth = 0

    # SMARTS patterns for Suzuki coupling components
    boronic_acid_pattern = Chem.MolFromSmarts("[B]([O])[O]")
    aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")

    def dfs_traverse(node, depth=0):
        nonlocal found_suzuki, suzuki_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if reactants match Suzuki coupling components
                has_boronic_acid = False
                has_aryl_halide = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(boronic_acid_pattern):
                            has_boronic_acid = True
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True

                if has_boronic_acid and has_aryl_halide:
                    found_suzuki = True
                    suzuki_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if Suzuki coupling is in the second half of the synthesis
    # (lower depth values are later in the synthesis)
    is_late_stage = suzuki_depth is not None and suzuki_depth < (max_depth / 2)
    result = found_suzuki and is_late_stage

    print(f"Late-stage Suzuki coupling strategy detected: {result}")
    print(f"  - Found Suzuki coupling: {found_suzuki}")
    print(f"  - Suzuki depth: {suzuki_depth}")
    print(f"  - Max depth: {max_depth}")
    print(f"  - Is late stage: {is_late_stage}")

    return result
