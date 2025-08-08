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
    This function detects a synthetic strategy involving an early Suzuki coupling
    for C-C bond formation between an aryl halide and a boronic acid/ester.
    """
    has_suzuki_coupling = False
    high_depth_reaction = False

    def dfs_traverse(node, depth=0):
        nonlocal has_suzuki_coupling, high_depth_reaction

        if node["type"] == "reaction":
            # Check if this is a high-depth (early in synthesis) reaction
            if depth >= 2:  # Depth 2 or higher is considered early in synthesis
                high_depth_reaction = True

                # Extract reactants and product
                rsmi = node["metadata"].get("rsmi", "")
                if rsmi:
                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]

                    reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
                    product = Chem.MolFromSmiles(product_part) if product_part else None

                    if product and all(r for r in reactants):
                        # Check for aryl halide
                        aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,Cl,I,F]")
                        # Check for boronic acid/ester
                        boronic_pattern = Chem.MolFromSmarts("[#6]B([#8])[#8]")

                        has_aryl_halide = any(
                            r.HasSubstructMatch(aryl_halide_pattern) for r in reactants
                        )
                        has_boronic = any(r.HasSubstructMatch(boronic_pattern) for r in reactants)

                        # Check if product has new C-C bond between aryl groups
                        if has_aryl_halide and has_boronic:
                            print("Detected potential Suzuki coupling at depth", depth)
                            has_suzuki_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    result = has_suzuki_coupling and high_depth_reaction
    print(f"Early Suzuki coupling strategy detected: {result}")
    return result
