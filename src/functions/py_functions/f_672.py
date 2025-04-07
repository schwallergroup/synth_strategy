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
    Detects if the synthetic route involves a nucleophilic aromatic substitution (SNAr)
    reaction to incorporate a morpholine group.
    """
    morpholine_snar_detected = False
    morpholine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#8][#6][#6]1")

    def dfs_traverse(node):
        nonlocal morpholine_snar_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for morpholine in reactants
                morpholine_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        morpholine_pattern
                    ):
                        morpholine_in_reactants = True
                        break

                # Check for halogenated aromatic in reactants
                halogenated_aromatic = False
                aromatic_halide_pattern = Chem.MolFromSmarts("[c]-[Cl,Br,I,F]")
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        aromatic_halide_pattern
                    ):
                        halogenated_aromatic = True
                        break

                # Check if product has morpholine attached to aromatic
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(morpholine_pattern):
                    if morpholine_in_reactants and halogenated_aromatic:
                        morpholine_snar_detected = True
                        print(f"Detected morpholine SNAr reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return morpholine_snar_detected
