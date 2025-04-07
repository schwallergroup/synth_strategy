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
    This function detects a synthetic strategy involving late-stage fragment coupling,
    specifically looking for complex fragments being joined in the final steps.
    """
    final_step_is_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_coupling

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a coupling reaction
            aryl_halide_pattern = Chem.MolFromSmarts("[c;R]-[#35,#53]")
            terminal_alkyne_pattern = Chem.MolFromSmarts("[C]#[CH]")
            diaryl_alkyne_pattern = Chem.MolFromSmarts("[c]-[C]#[C]-[c]")

            # Check reactants for aryl halide and alkyne
            has_aryl_halide = False
            has_alkyne = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True
                    if mol.HasSubstructMatch(terminal_alkyne_pattern):
                        has_alkyne = True

            # Check product for diaryl alkyne
            product_mol = Chem.MolFromSmiles(product)
            has_coupled_product = False
            if product_mol and product_mol.HasSubstructMatch(diaryl_alkyne_pattern):
                has_coupled_product = True

            is_coupling = has_aryl_halide and has_alkyne and has_coupled_product

            # Check if this is a late-stage coupling (depth 0 or 1)
            if depth <= 1 and is_coupling:
                final_step_is_coupling = True
                print(f"Found late-stage coupling reaction at depth {depth}: {rsmi}")

            # Check if fragments are complex (contain multiple rings or functional groups)
            if is_coupling:
                complex_fragments = 0
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        ring_count = Chem.GetSSSR(mol)
                        if len(ring_count) >= 1:
                            complex_fragments += 1

                if complex_fragments >= 2:
                    print(f"Found coupling of complex fragments: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage fragment coupling strategy detected: {final_step_is_coupling}")
    return final_step_is_coupling
