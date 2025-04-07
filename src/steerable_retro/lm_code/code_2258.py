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
    Detects if the synthesis involves functional group interconversion (particularly
    amine to halide) to prepare a substrate for cross-coupling.
    """
    # Track if we found relevant transformations
    found_halogenation = False
    found_coupling = False
    halogenation_depth = None
    coupling_depth = None

    def is_amine_to_halide(reaction_smiles):
        """Check if a reaction converts an amine to a halide"""
        # Split reaction into reactants and products
        parts = reaction_smiles.split(">")
        if len(parts) != 3:
            return False

        reactants = parts[0].split(".")
        product = parts[2]

        # Check for amine in reactants
        has_amine = False
        for reactant in reactants:
            try:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[N]")):
                    has_amine = True
                    break
            except:
                continue

        # Check for aryl halide in product
        try:
            prod_mol = Chem.MolFromSmiles(product)
            has_aryl_halide = prod_mol and prod_mol.HasSubstructMatch(
                Chem.MolFromSmarts("[c]-[#53,#35,#17]")
            )

            # If reactant has amine and product has halide, likely an amine to halide conversion
            if has_amine and has_aryl_halide:
                return True
        except:
            pass

        return False

    def is_cross_coupling(reaction_smiles):
        """Check if a reaction is a cross-coupling reaction"""
        # Patterns for cross-coupling
        aryl_halide_pattern = Chem.MolFromSmarts("[c]-[#53,#35,#17]")  # Aryl-I, Br, Cl

        # Split reaction into reactants and products
        parts = reaction_smiles.split(">")
        if len(parts) != 3:
            return False

        reactants = parts[0].split(".")
        product = parts[2]

        # Check if reactants contain aryl halide
        has_aryl_halide = False
        for reactant in reactants:
            try:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                    has_aryl_halide = True
                    break
            except:
                continue

        # Check if product has a new C-C bond
        try:
            prod_mol = Chem.MolFromSmiles(product)
            if prod_mol and has_aryl_halide:
                # If reactant has aryl halide and product has new bonds, likely a cross-coupling
                return True
        except:
            pass

        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_halogenation, found_coupling, halogenation_depth, coupling_depth

        if node["type"] == "reaction":
            # Get reaction SMILES
            if "rsmi" in node.get("metadata", {}):
                reaction_smiles = node["metadata"]["rsmi"]

                # Check if this is an amine to halide conversion
                if is_amine_to_halide(reaction_smiles):
                    found_halogenation = True
                    halogenation_depth = depth
                    print(f"Found amine to halide conversion at depth {depth}")

                # Check if this is a cross-coupling
                if is_cross_coupling(reaction_smiles):
                    found_coupling = True
                    coupling_depth = depth
                    print(f"Found cross-coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found both transformations in the correct order
    if (
        found_halogenation
        and found_coupling
        and halogenation_depth is not None
        and coupling_depth is not None
    ):
        # Halogenation should occur before coupling in the synthesis (higher depth in retrosynthesis)
        if halogenation_depth > coupling_depth:
            print("Detected functional group interconversion for coupling strategy")
            return True

    return False
