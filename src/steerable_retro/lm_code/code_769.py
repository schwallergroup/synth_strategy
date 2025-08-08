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
    This function detects Stille coupling reactions (organotin reagent + aryl halide → biaryl).
    """
    stille_detected = False

    def dfs_traverse(node):
        nonlocal stille_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for organotin reagent in reactants
            tin_pattern = Chem.MolFromSmarts("[#50]")  # Tin atom

            # Check for aryl halide in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("c[Br,I,Cl]")

            # Check for biaryl in product
            biaryl_pattern = Chem.MolFromSmarts("c-c")  # Simplified biaryl pattern

            tin_present = False
            aryl_halide_present = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(tin_pattern):
                        tin_present = True
                        print(f"Found tin reagent: {reactant}")
                    if mol.HasSubstructMatch(aryl_halide_pattern):
                        aryl_halide_present = True
                        print(f"Found aryl halide: {reactant}")

            product_mol = Chem.MolFromSmiles(product)
            biaryl_present = False
            if product_mol and product_mol.HasSubstructMatch(biaryl_pattern):
                biaryl_present = True
                print(f"Found biaryl in product: {product}")

            if tin_present and aryl_halide_present and biaryl_present:
                stille_detected = True
                print("Stille coupling detected!")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return stille_detected
