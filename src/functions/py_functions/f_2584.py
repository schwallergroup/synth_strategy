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
    Detects if the synthesis uses a late-stage thiazole formation strategy,
    where a thiazole ring is formed in the final step from a thiourea intermediate.
    """
    # Track if we found thiourea formation and thiazole formation
    thiourea_formation_found = False
    thiazole_formation_found = False
    thiazole_formation_depth = float("inf")
    thiourea_formation_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal thiourea_formation_found, thiazole_formation_found
        nonlocal thiazole_formation_depth, thiourea_formation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for thiourea formation
            thiourea_pattern = Chem.MolFromSmarts("[#6](=[#16])[#7]")
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and product_mol.HasSubstructMatch(thiourea_pattern):
                # Check if any reactant doesn't have thiourea
                reactant_has_thiourea = False
                for reactant_smiles in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        thiourea_pattern
                    ):
                        reactant_has_thiourea = True
                        break

                if not reactant_has_thiourea:
                    thiourea_formation_found = True
                    thiourea_formation_depth = depth
                    print(f"Found thiourea formation at depth {depth}")

            # Check for thiazole formation
            thiazole_pattern = Chem.MolFromSmarts("[s]1[c][n][c][c]1")
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and product_mol.HasSubstructMatch(thiazole_pattern):
                # Check if any reactant doesn't have thiazole
                reactant_has_thiazole = False
                for reactant_smiles in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        thiazole_pattern
                    ):
                        reactant_has_thiazole = True
                        break

                if not reactant_has_thiazole:
                    thiazole_formation_found = True
                    thiazole_formation_depth = depth
                    print(f"Found thiazole formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have both thiourea formation and thiazole formation,
    # and thiazole formation happens at a lower depth (later in synthesis)
    # than thiourea formation
    if (
        thiourea_formation_found
        and thiazole_formation_found
        and thiazole_formation_depth < thiourea_formation_depth
    ):
        print("Detected late-stage thiazole formation strategy")
        return True
    return False
