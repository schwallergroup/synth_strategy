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
    This function detects if the synthetic route includes nucleophilic aromatic substitution
    with a piperazine displacing fluorine on a nitro-activated aromatic ring.
    """
    has_snar = False

    def dfs_traverse(node):
        nonlocal has_snar

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl fluoride with nitro group in reactants
            aryl_f_nitro_pattern = Chem.MolFromSmarts(
                "[F]-[c]1[c]([N+](=O)[O-])[c][c][c][c]1"
            )
            piperazine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6][#6]1")

            has_aryl_f_nitro = False
            has_piperazine = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    if reactant_mol.HasSubstructMatch(aryl_f_nitro_pattern):
                        has_aryl_f_nitro = True
                    if reactant_mol.HasSubstructMatch(piperazine_pattern):
                        has_piperazine = True

            # Check for C-N bond formation
            if has_aryl_f_nitro and has_piperazine:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Check if product has piperazine attached to aromatic ring where F was
                    c_n_pattern = Chem.MolFromSmarts("[c]-[#7]1[#6][#6][#7][#6][#6]1")
                    if product_mol.HasSubstructMatch(c_n_pattern):
                        has_snar = True
                        print(f"Nucleophilic aromatic substitution detected: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Nucleophilic aromatic substitution detected: {has_snar}")
    return has_snar
