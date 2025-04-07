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
    Detects convergent synthesis with late-stage SNAr coupling of a chloropyrimidine and a phenol.
    """
    # Track if we found the SNAr coupling
    found_snar_coupling = False

    def dfs_traverse(node):
        nonlocal found_snar_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is the final coupling reaction (depth 0)
            if len(node.get("children", [])) == 2:  # Has two reactant children
                # Check for chloropyrimidine pattern in one reactant
                chloropyrimidine_pattern = Chem.MolFromSmarts("[c]1[n][c][c][n][c]1[Cl]")
                # Check for phenol pattern in the other reactant
                phenol_pattern = Chem.MolFromSmarts("[OH][c]")

                # Check for ether bond formation in product
                ether_pyrimidine_pattern = Chem.MolFromSmarts("[c]1[n][c][c][n][c]1[O][c]")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and all(r for r in reactant_mols):
                    has_chloropyrimidine = any(
                        m.HasSubstructMatch(chloropyrimidine_pattern) for m in reactant_mols if m
                    )
                    has_phenol = any(
                        m.HasSubstructMatch(phenol_pattern) for m in reactant_mols if m
                    )
                    forms_ether = product_mol.HasSubstructMatch(ether_pyrimidine_pattern)

                    if has_chloropyrimidine and has_phenol and forms_ether:
                        print("Found SNAr coupling of chloropyrimidine and phenol")
                        found_snar_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_snar_coupling
