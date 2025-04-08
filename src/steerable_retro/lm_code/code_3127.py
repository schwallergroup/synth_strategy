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
    This function detects a strategy involving aromatic nucleophilic substitution
    where a piperidine connects to a halogenated pyridine.
    """
    # Track if we find the pattern
    found_pattern = False

    # SMARTS patterns
    chloropyridine = Chem.MolFromSmarts("[n]1[cH][c]([Cl])[c][cH][c]1")
    n_methylpiperidine = Chem.MolFromSmarts("[CH3][N]1[CH2][CH2][CH][CH2][CH2]1")
    piperidine_pyridine_connection = Chem.MolFromSmarts(
        "[CH3][N]1[CH2][CH2][CH]([NH][c]2[cH][c][n][cH][c]2)[CH2][CH2]1"
    )

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for reactants with chloropyridine and piperidine
            has_chloropyridine = False
            has_piperidine = False

            for reactant_smiles in reactants_smiles:
                try:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(chloropyridine):
                            has_chloropyridine = True
                        if reactant_mol.HasSubstructMatch(n_methylpiperidine):
                            has_piperidine = True
                except:
                    continue

            # Check if product has the connection
            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(piperidine_pyridine_connection):
                    if has_chloropyridine and has_piperidine:
                        print(f"Found aromatic nucleophilic substitution at depth {depth}")
                        found_pattern = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_pattern
