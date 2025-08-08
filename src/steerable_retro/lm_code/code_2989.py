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
    This function detects if benzyl ether formation is used as a protection strategy.
    """
    has_benzyl_protection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_benzyl_protection

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains benzyl ether
                product_mol = Chem.MolFromSmiles(product)
                benzyl_ether_pattern = Chem.MolFromSmarts("[c][CH2][O][c]")

                # Check if reactants contain phenol and benzyl halide
                has_phenol = False
                has_benzyl_halide = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        phenol_pattern = Chem.MolFromSmarts("[c][OH]")
                        benzyl_halide_pattern = Chem.MolFromSmarts("[c][CH2][Br,Cl,I,F]")

                        if reactant_mol.HasSubstructMatch(phenol_pattern):
                            has_phenol = True
                        if reactant_mol.HasSubstructMatch(benzyl_halide_pattern):
                            has_benzyl_halide = True

                # If product has benzyl ether and reactants have phenol and benzyl halide
                if (
                    product_mol
                    and product_mol.HasSubstructMatch(benzyl_ether_pattern)
                    and has_phenol
                    and has_benzyl_halide
                ):
                    has_benzyl_protection = True
                    print(f"Detected benzyl ether protection at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return has_benzyl_protection
