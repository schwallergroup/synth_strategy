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
    Detects a borylation followed by Suzuki coupling sequence.
    """
    borylation_detected = False
    suzuki_detected = False
    borylation_depth = -1
    suzuki_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal borylation_detected, suzuki_detected, borylation_depth, suzuki_depth

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                # Check for borylation (aryl-Br to aryl-B(OR)2)
                if product_mol is not None:
                    boronate_pattern = Chem.MolFromSmarts("[c]-[B;X3]")
                    if product_mol.HasSubstructMatch(boronate_pattern):
                        halide_pattern = Chem.MolFromSmarts("[c]-[Br,I,Cl]")
                        if any(
                            mol is not None and mol.HasSubstructMatch(halide_pattern)
                            for mol in reactant_mols
                        ):
                            borylation_detected = True
                            borylation_depth = depth
                            print(f"Detected borylation at depth {depth}")

                # Check for Suzuki coupling
                if product_mol is not None:
                    boronate_pattern = Chem.MolFromSmarts("[c]-[B;X3]")
                    halide_pattern = Chem.MolFromSmarts("[c]-[Br,I,Cl]")

                    has_boronate = any(
                        mol is not None and mol.HasSubstructMatch(boronate_pattern)
                        for mol in reactant_mols
                    )
                    has_halide = any(
                        mol is not None and mol.HasSubstructMatch(halide_pattern)
                        for mol in reactant_mols
                    )

                    if has_boronate and has_halide:
                        suzuki_detected = True
                        suzuki_depth = depth
                        print(f"Detected Suzuki coupling at depth {depth}")

            except Exception as e:
                print(f"Error in borylation-Suzuki detection: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if borylation occurs before Suzuki coupling
    sequence_detected = (
        borylation_detected and suzuki_detected and borylation_depth > suzuki_depth
    )  # Remember: higher depth = earlier in synthesis

    return sequence_detected
