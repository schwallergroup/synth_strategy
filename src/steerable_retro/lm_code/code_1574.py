#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    Detects if the synthesis route involves SNAr with morpholine as a nucleophile.
    """
    has_snar_morpholine = False

    def dfs_traverse(node, depth=0):
        nonlocal has_snar_morpholine

        if node["type"] == "reaction":
            # Get reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for morpholine in reactants
            morpholine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#8][#6][#6]1")

            # Check for aryl halide in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("c-[#17,#35,#53]")

            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(reactant_mols):
                has_morpholine = any(
                    [mol and mol.HasSubstructMatch(morpholine_pattern) for mol in reactant_mols]
                )
                has_aryl_halide = any(
                    [mol and mol.HasSubstructMatch(aryl_halide_pattern) for mol in reactant_mols]
                )

                # Check if product has morpholine connected to aryl
                aryl_morpholine_pattern = Chem.MolFromSmarts("c-[#7]1[#6][#6][#8][#6][#6]1")
                has_aryl_morpholine = product_mol.HasSubstructMatch(aryl_morpholine_pattern)

                if has_morpholine and has_aryl_halide and has_aryl_morpholine:
                    has_snar_morpholine = True
                    print(f"Detected SNAr with morpholine at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_snar_morpholine
