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
    This function detects if the synthesis includes a convergent approach where two complex heterocyclic fragments are combined.
    """
    convergent_assembly_detected = False

    def dfs_traverse(node):
        nonlocal convergent_assembly_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Need at least 2 reactants for convergent assembly
            if len(reactants) >= 2:
                # Check if both reactants contain heterocycles
                heterocycle_patterns = [
                    Chem.MolFromSmarts("[#6]1:[#7]:[#6]:[#6]:[#6]:[#6]:1"),  # pyridine
                    Chem.MolFromSmarts("[#6]1:[#7]:[#6]:[#7]:[#6]:1"),  # pyrimidine
                    Chem.MolFromSmarts("[#6]1:[#7]:[#6]:[#6]:[#7]:1"),  # pyrazine
                    Chem.MolFromSmarts("[#6]1:[#7]:[#7]:[#6]:[#6]:1"),  # pyridazine
                    Chem.MolFromSmarts(
                        "[#6]1:[#7]:[#6]:[#6]:[#6]:[#7]:1"
                    ),  # quinoxaline
                ]

                complex_heterocycles = []
                for i, reactant in enumerate(reactants):
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    # Count atoms to determine complexity
                    atom_count = reactant_mol.GetNumAtoms()

                    # Check for heterocycle presence
                    for pattern in heterocycle_patterns:
                        if reactant_mol.HasSubstructMatch(pattern) and atom_count > 8:
                            complex_heterocycles.append(i)
                            break

                # If we have at least 2 complex heterocyclic reactants
                if len(complex_heterocycles) >= 2:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # Check if product is more complex than individual reactants
                        product_atom_count = product_mol.GetNumAtoms()
                        if product_atom_count > max(
                            [
                                Chem.MolFromSmiles(reactants[i]).GetNumAtoms()
                                for i in complex_heterocycles
                            ]
                        ):
                            convergent_assembly_detected = True
                            print("Detected convergent heterocycle assembly")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return convergent_assembly_detected
