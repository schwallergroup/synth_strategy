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
    Detects if the synthesis involves an esterification followed by
    conversion to a hydrazide.
    """
    has_ester_formation = False
    has_hydrazide_formation = False
    ester_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal has_ester_formation, has_hydrazide_formation, ester_depth

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for esterification
            has_carboxylic_acid = False
            has_alcohol = False

            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
                    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H1]")

                    if mol.HasSubstructMatch(acid_pattern):
                        has_carboxylic_acid = True
                    if mol.HasSubstructMatch(alcohol_pattern):
                        has_alcohol = True

            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol:
                # Check for ester formation
                if has_carboxylic_acid and has_alcohol:
                    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][CX4]")
                    if product_mol.HasSubstructMatch(ester_pattern):
                        has_ester_formation = True
                        ester_depth = depth
                        print(f"Ester formation detected at depth {depth}")

                # Check for hydrazide formation
                hydrazine_pattern = Chem.MolFromSmarts("[NX3][NX3]")
                for reactant in reactants_smiles:
                    if Chem.MolFromSmiles(reactant) and Chem.MolFromSmiles(
                        reactant
                    ).HasSubstructMatch(hydrazine_pattern):
                        hydrazide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3][NX3]")
                        if product_mol.HasSubstructMatch(hydrazide_pattern):
                            has_hydrazide_formation = True
                            hydrazide_depth = depth
                            print(f"Hydrazide formation detected at depth {depth}")

                            # Check if this follows esterification
                            if has_ester_formation and ester_depth > hydrazide_depth:
                                print("Esterification-hydrazide sequence detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    # Return true if both reactions occurred and in the correct sequence
    return has_ester_formation and has_hydrazide_formation and ester_depth > 0
