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
    This function detects biaryl formation via coupling reactions.
    """
    biaryl_formation_found = False

    def dfs_traverse(node):
        nonlocal biaryl_formation_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for biaryl formation
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # SMARTS pattern for biaryl
                    biaryl_pattern = Chem.MolFromSmarts("c1ccccc1-c1ccccc1")
                    if product_mol.HasSubstructMatch(biaryl_pattern):
                        # Check if reactants don't have biaryl
                        reactants_have_biaryl = False
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                biaryl_pattern
                            ):
                                reactants_have_biaryl = True
                                break

                        if not reactants_have_biaryl:
                            # Check for boronic acid or halide in reactants (typical for Suzuki coupling)
                            boronic_acid_pattern = Chem.MolFromSmarts(
                                "[#6]-[#5](-[#8])-[#8]"
                            )
                            halide_pattern = Chem.MolFromSmarts("c-[#17,#35,#53]")

                            has_boronic_acid = False
                            has_halide = False

                            for reactant in reactants:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if not reactant_mol:
                                    continue

                                if reactant_mol.HasSubstructMatch(boronic_acid_pattern):
                                    has_boronic_acid = True
                                if reactant_mol.HasSubstructMatch(halide_pattern):
                                    has_halide = True

                            if has_boronic_acid and has_halide:
                                print("Biaryl formation via coupling detected")
                                biaryl_formation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return biaryl_formation_found
