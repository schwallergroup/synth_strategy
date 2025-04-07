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
    Detects if the synthetic route involves a borylation reaction followed by a coupling reaction.
    """
    borylation_depths = []
    coupling_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for borylation: Ar-X → Ar-B(OR)2
                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,I,Cl]")
                boronic_pattern = Chem.MolFromSmarts("[c]-[B]([O])[O]")
                boronic_ester_pattern = Chem.MolFromSmarts("[c]-[B]1[O][C][C][O]1")

                # Check for coupling: formation of biaryl
                biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

                # Check for borylation
                has_aryl_halide_reactant = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide_reactant = True
                        break

                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol:
                    # Check for borylation
                    if has_aryl_halide_reactant and (
                        prod_mol.HasSubstructMatch(boronic_pattern)
                        or prod_mol.HasSubstructMatch(boronic_ester_pattern)
                    ):
                        print(f"Detected borylation reaction at depth {depth}")
                        borylation_depths.append(depth)

                    # Check for coupling
                    if prod_mol.HasSubstructMatch(biaryl_pattern):
                        # Check if one reactant has boronic group
                        has_boronic_reactant = False
                        for reactant in reactants:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and (
                                mol.HasSubstructMatch(boronic_pattern)
                                or mol.HasSubstructMatch(boronic_ester_pattern)
                            ):
                                has_boronic_reactant = True
                                break

                        if has_boronic_reactant:
                            print(f"Detected coupling reaction at depth {depth}")
                            coupling_depths.append(depth)

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if there's a borylation followed by a coupling
    for b_depth in borylation_depths:
        for c_depth in coupling_depths:
            if c_depth < b_depth:  # Remember: lower depth = later in synthesis
                print(
                    f"Detected borylation at depth {b_depth} followed by coupling at depth {c_depth}"
                )
                return True

    return False
