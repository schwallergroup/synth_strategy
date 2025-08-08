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
    This function detects metal-catalyzed coupling reactions (like Stille coupling)
    for C-C bond formation, typically involving an organometallic reagent and aryl halide.
    """
    has_metal_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_metal_coupling

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for organometallic reagent (tin compound for Stille)
            has_organometallic = False
            has_aryl_halide = False

            for reactant in reactants_smiles:
                if reactant and Chem.MolFromSmiles(reactant):
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    # Check for tin (Stille coupling)
                    if "[Sn]" in reactant:
                        has_organometallic = True
                        print(f"Found organometallic reagent in reaction at depth {depth}")

                    # Check for aryl halide
                    aryl_halide_pattern = Chem.MolFromSmarts("c[Br,I,Cl]")
                    if reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True
                        print(f"Found aryl halide in reaction at depth {depth}")

            # If both organometallic and aryl halide are present, likely a coupling reaction
            if has_organometallic and has_aryl_halide:
                has_metal_coupling = True
                print(f"Detected metal-catalyzed coupling reaction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_metal_coupling
