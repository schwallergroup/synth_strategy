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
    This function detects a convergent synthesis strategy that includes a Suzuki coupling
    to form a biaryl motif.
    """
    # Track if we found the strategy components
    found_suzuki = False
    suzuki_depth = None
    found_convergent = False

    # SMARTS patterns
    aryl_halide_pattern = Chem.MolFromSmarts("c1ccccc1[I,Br,Cl]")
    boronic_acid_pattern = Chem.MolFromSmarts("c1ccccc1B(O)O")

    def dfs_traverse(node, depth=0):
        nonlocal found_suzuki, suzuki_depth, found_convergent

        if node["type"] == "reaction":
            # Check if this is a Suzuki coupling
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if we have exactly 2 reactants (typical for Suzuki)
                if len(reactants) == 2:
                    try:
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                        product_mol = Chem.MolFromSmiles(product)

                        # Check for aryl halide and boronic acid patterns
                        has_aryl_halide = any(
                            mol is not None and mol.HasSubstructMatch(aryl_halide_pattern)
                            for mol in reactant_mols
                        )
                        has_boronic_acid = any(
                            mol is not None and mol.HasSubstructMatch(boronic_acid_pattern)
                            for mol in reactant_mols
                        )

                        if has_aryl_halide and has_boronic_acid and product_mol is not None:
                            print(f"Found Suzuki coupling at depth {depth}")
                            found_suzuki = True
                            suzuki_depth = depth
                    except:
                        print("Error processing SMILES in reaction")

            # Check if this is a convergent step (combining two complex fragments)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If we have multiple reactants and we're in the middle of the synthesis
                if len(reactants) >= 2 and 0 < depth < 3:
                    found_convergent = True
                    print(f"Found convergent step at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we found both Suzuki coupling and convergent synthesis
    strategy_present = found_suzuki and found_convergent
    print(f"Convergent synthesis with Suzuki coupling: {strategy_present}")
    return strategy_present
