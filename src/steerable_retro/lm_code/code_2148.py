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
    This function detects the introduction of a morpholine group via nucleophilic substitution.
    """
    morpholine_introduced = False

    def dfs_traverse(node):
        nonlocal morpholine_introduced

        if node["type"] == "reaction":
            # Check if this is a reaction node with metadata
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alkyl halide and morpholine in reactants
                alkyl_halide_pattern = Chem.MolFromSmarts("[#6][Cl,Br,I,F]")
                morpholine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#8][#6][#6]1")

                has_alkyl_halide = False
                has_morpholine = False

                for reactant in reactants:
                    if reactant:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                if mol.HasSubstructMatch(alkyl_halide_pattern):
                                    has_alkyl_halide = True
                                if mol.HasSubstructMatch(morpholine_pattern):
                                    has_morpholine = True
                        except:
                            continue

                # Check for alkylated morpholine in product
                if has_alkyl_halide and has_morpholine:
                    alkylated_morpholine_pattern = Chem.MolFromSmarts(
                        "[#6][#7]1[#6][#6][#8][#6][#6]1"
                    )
                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(
                            alkylated_morpholine_pattern
                        ):
                            print("Detected morpholine introduction")
                            morpholine_introduced = True
                    except:
                        pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return morpholine_introduced
