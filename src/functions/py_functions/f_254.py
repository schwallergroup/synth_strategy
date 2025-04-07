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
    This function detects if the synthesis involves sulfonamide formation
    followed by N-alkylation of the sulfonamide
    """
    # Track sulfonamide formation and alkylation reactions and their depths
    sulfonamide_formation_depth = -1
    sulfonamide_alkylation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_formation_depth, sulfonamide_alkylation_depth

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for sulfonamide formation
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts("[NH]S(=O)(=O)[#6]")
            ):

                # Check if reactants include sulfonyl chloride and amine
                has_sulfonyl_chloride = False
                has_amine = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[#6]S(=O)(=O)Cl")
                        ):
                            has_sulfonyl_chloride = True
                        if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
                            has_amine = True

                if has_sulfonyl_chloride or has_amine:
                    sulfonamide_formation_depth = depth
                    print(f"Found sulfonamide formation at depth {depth}")

            # Check for N-alkylation of sulfonamide
            if (
                product_mol
                and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[#6][N]S(=O)(=O)[#6]")
                )
                and not product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[NH]S(=O)(=O)[#6]")
                )
            ):

                # Check if reactants include secondary sulfonamide
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[NH]S(=O)(=O)[#6]")
                    ):
                        sulfonamide_alkylation_depth = depth
                        print(f"Found sulfonamide N-alkylation at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if both transformations were found in the correct order
    # Remember that lower depth values are later in the synthesis
    correct_sequence = (
        sulfonamide_formation_depth != -1
        and sulfonamide_alkylation_depth != -1
        and sulfonamide_formation_depth > sulfonamide_alkylation_depth
    )

    print(f"Sulfonamide formation and alkylation strategy detected: {correct_sequence}")
    return correct_sequence
