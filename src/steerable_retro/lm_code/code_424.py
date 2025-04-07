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
    This function detects if the synthesis involves a sequence of transformations
    from nitro group to sulfonamide (nitro → amine → sulfonamide).
    """
    nitro_reduction_depth = None
    sulfonamide_formation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_depth, sulfonamide_formation_depth

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro reduction
            nitro_in_reactants = False
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[N+](=[O])[O-]")
                ):
                    nitro_in_reactants = True
                    break

            product_mol = Chem.MolFromSmiles(product_smiles)
            if (
                nitro_in_reactants
                and product_mol
                and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]"))
            ):
                nitro_reduction_depth = depth
                print(f"Nitro reduction detected at depth {depth}")

            # Check for sulfonamide formation
            amine_in_reactants = False
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
                    amine_in_reactants = True
                    break

            if (
                amine_in_reactants
                and product_mol
                and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[N][S](=O)(=O)[C]"))
            ):
                sulfonamide_formation_depth = depth
                print(f"Sulfonamide formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if both transformations were found and in the correct order
    has_sequence = (
        nitro_reduction_depth is not None
        and sulfonamide_formation_depth is not None
        and nitro_reduction_depth > sulfonamide_formation_depth
    )

    print(f"Nitro to sulfonamide sequence: {has_sequence}")
    print(
        f"Nitro reduction depth: {nitro_reduction_depth}, Sulfonamide formation depth: {sulfonamide_formation_depth}"
    )
    return has_sequence
