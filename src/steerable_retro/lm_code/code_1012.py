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
    Detects a convergent synthesis strategy using alkyne coupling chemistry.
    Looks for:
    1. Presence of alkyne coupling reaction (C(sp)-C(sp²))
    2. Multiple fragments being combined
    3. TMS-alkyne deprotection
    """
    # Track if we found the key features
    found_alkyne_coupling = False
    found_tms_deprotection = False
    fragment_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal found_alkyne_coupling, found_tms_deprotection, fragment_count

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Count fragments being combined
                if len(reactants) > 1:
                    fragment_count += 1
                    print(
                        f"Found fragment combination at depth {depth}: {len(reactants)} fragments"
                    )

                # Check for alkyne coupling (aryl halide + alkyne)
                aryl_halide_pattern = Chem.MolFromSmarts("c-[Br,I,Cl]")
                terminal_alkyne_pattern = Chem.MolFromSmarts("[C]#[C;H]")
                internal_alkyne_pattern = Chem.MolFromSmarts("[C]#[C;!H]")

                product_mol = Chem.MolFromSmiles(product)

                # Check if product has internal alkyne connected to aryl
                if product_mol and product_mol.HasSubstructMatch(internal_alkyne_pattern):
                    # Check if reactants include aryl halide and terminal alkyne
                    has_aryl_halide = False
                    has_terminal_alkyne = False

                    for reactant in reactants:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol:
                            if r_mol.HasSubstructMatch(aryl_halide_pattern):
                                has_aryl_halide = True
                            if r_mol.HasSubstructMatch(terminal_alkyne_pattern):
                                has_terminal_alkyne = True

                    if has_aryl_halide and has_terminal_alkyne:
                        found_alkyne_coupling = True
                        print(f"Found alkyne coupling at depth {depth}")

                # Check for TMS-alkyne deprotection
                tms_alkyne_pattern = Chem.MolFromSmarts("[C]-[Si](-[C])(-[C])-[C]")

                for reactant in reactants:
                    r_mol = Chem.MolFromSmiles(reactant)
                    if r_mol and r_mol.HasSubstructMatch(tms_alkyne_pattern):
                        if product_mol and product_mol.HasSubstructMatch(terminal_alkyne_pattern):
                            found_tms_deprotection = True
                            print(f"Found TMS-alkyne deprotection at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if key features were found
    return found_alkyne_coupling and found_tms_deprotection and fragment_count >= 1
