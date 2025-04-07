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
    This function detects a synthetic strategy involving orthogonal modifications of
    different functional groups (side chain vs. aromatic substituents).
    """
    side_chain_modifications = 0
    aromatic_modifications = 0

    def dfs_traverse(node):
        nonlocal side_chain_modifications, aromatic_modifications

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for side chain modifications (alcohol/aldehyde/acid)
            aldehyde_pattern = Chem.MolFromSmarts("[C;H1](=O)")
            alcohol_pattern = Chem.MolFromSmarts("[C;H1,H2]([OH])")
            carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")

            # Check for aromatic substituent modifications
            aryl_halide_pattern = Chem.MolFromSmarts("[c]-[#9,#17,#35,#53]")
            aniline_pattern = Chem.MolFromSmarts("[c]-[N;H2]")
            nitro_pattern = Chem.MolFromSmarts("[c]-[N+](=O)[O-]")

            r_mol = Chem.MolFromSmiles(reactants[0])
            p_mol = Chem.MolFromSmiles(product)

            if r_mol and p_mol:
                # Check side chain modifications
                r_has_aldehyde = r_mol.HasSubstructMatch(aldehyde_pattern)
                r_has_alcohol = r_mol.HasSubstructMatch(alcohol_pattern)
                r_has_acid = r_mol.HasSubstructMatch(carboxylic_acid_pattern)

                p_has_aldehyde = p_mol.HasSubstructMatch(aldehyde_pattern)
                p_has_alcohol = p_mol.HasSubstructMatch(alcohol_pattern)
                p_has_acid = p_mol.HasSubstructMatch(carboxylic_acid_pattern)

                if (
                    r_has_aldehyde != p_has_aldehyde
                    or r_has_alcohol != p_has_alcohol
                    or r_has_acid != p_has_acid
                ):
                    side_chain_modifications += 1
                    print(
                        f"Detected side chain modification: {reactants[0]} -> {product}"
                    )

                # Check aromatic modifications
                r_has_halide = r_mol.HasSubstructMatch(aryl_halide_pattern)
                r_has_amine = r_mol.HasSubstructMatch(aniline_pattern)
                r_has_nitro = r_mol.HasSubstructMatch(nitro_pattern)

                p_has_halide = p_mol.HasSubstructMatch(aryl_halide_pattern)
                p_has_amine = p_mol.HasSubstructMatch(aniline_pattern)
                p_has_nitro = p_mol.HasSubstructMatch(nitro_pattern)

                if (
                    r_has_halide != p_has_halide
                    or r_has_amine != p_has_amine
                    or r_has_nitro != p_has_nitro
                ):
                    aromatic_modifications += 1
                    print(
                        f"Detected aromatic modification: {reactants[0]} -> {product}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both side chain and aromatic modifications are present
    return side_chain_modifications > 0 and aromatic_modifications > 0
