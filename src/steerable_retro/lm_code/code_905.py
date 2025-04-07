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
    This function detects if the synthesis involves sequential aromatic coupling reactions.
    """
    coupling_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal coupling_depths

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for boronic acid/ester pattern in reactants
            boronic_pattern = Chem.MolFromSmarts("[c][B]([O])[O]")
            boronic_ester_pattern = Chem.MolFromSmarts("[c][B]1[O][C]([C])([C])[C]([C])([C])[O]1")

            # Check for aryl halide pattern in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")

            # Check for biaryl formation in product
            biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

            has_boronic = False
            has_aryl_halide = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(boronic_pattern) or mol.HasSubstructMatch(
                            boronic_ester_pattern
                        ):
                            has_boronic = True
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                has_biaryl = product_mol and product_mol.HasSubstructMatch(biaryl_pattern)
            except:
                has_biaryl = False

            if has_boronic and has_aryl_halide and has_biaryl:
                print(f"Found aromatic coupling at depth {depth}: {rsmi}")
                coupling_depths.append(depth)

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if there are at least 2 coupling reactions with different depths
    sequential = len(set(coupling_depths)) >= 2
    print(f"Sequential aromatic couplings found: {sequential} at depths {coupling_depths}")
    return sequential
