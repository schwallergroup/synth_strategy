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
    Detects if the synthetic route involves a late-stage Suzuki coupling (at depth 0 or 1).
    """
    late_stage_suzuki = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_suzuki

        if node["type"] == "reaction" and depth <= 1:  # Only check reactions at depth 0 or 1
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid/ester in one reactant
                boronic_pattern = Chem.MolFromSmarts("[c]-[B]([O])[O]")
                boronic_ester_pattern = Chem.MolFromSmarts("[c]-[B]1[O][C][C][O]1")

                # Check for aryl halide in another reactant
                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,I,Cl]")

                has_boronic = False
                has_aryl_halide = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(boronic_pattern) or mol.HasSubstructMatch(
                            boronic_ester_pattern
                        ):
                            has_boronic = True
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True

                # Check if product has biaryl system
                if has_boronic and has_aryl_halide:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol:
                        # Simple check for biaryl formation
                        biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")
                        if prod_mol.HasSubstructMatch(biaryl_pattern):
                            print(f"Detected late-stage Suzuki coupling at depth {depth}")
                            late_stage_suzuki = True

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_suzuki
