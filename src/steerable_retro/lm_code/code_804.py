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
    This function detects if the synthesis uses a sequential strategy of
    halogenation followed by coupling reaction.
    """
    halogenation_depths = []
    coupling_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for halogenation (specifically iodination)
            iodine_pattern = Chem.MolFromSmarts("[#6][#53]")
            product_mol = Chem.MolFromSmiles(product)

            if product_mol and product_mol.HasSubstructMatch(iodine_pattern):
                # Check if iodine wasn't in reactants
                iodine_in_reactants = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(iodine_pattern):
                        iodine_in_reactants = True
                        break

                if not iodine_in_reactants:
                    print(f"Halogenation (iodination) detected at depth {depth}")
                    halogenation_depths.append(depth)

            # Check for coupling reaction (simplified)
            boronic_acid_pattern = Chem.MolFromSmarts("[#6]B(O)O")
            aryl_halide_pattern = Chem.MolFromSmarts("[#6]~[#53,#35,#17]")

            boronic_acid_present = False
            aryl_halide_present = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(boronic_acid_pattern):
                        boronic_acid_present = True
                    if mol.HasSubstructMatch(aryl_halide_pattern):
                        aryl_halide_present = True

            if boronic_acid_present and aryl_halide_present:
                print(f"Coupling reaction detected at depth {depth}")
                coupling_depths.append(depth)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if halogenation is followed by coupling
    for h_depth in halogenation_depths:
        for c_depth in coupling_depths:
            # In retrosynthetic direction, coupling should be at lower depth than halogenation
            if c_depth < h_depth:
                print(
                    f"Sequential halogenation-coupling strategy detected: halogenation at depth {h_depth}, coupling at depth {c_depth}"
                )
                return True

    return False
