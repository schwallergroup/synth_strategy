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
    Detects if the synthesis involves sequential functionalization of aryl halides
    (e.g., nitration followed by cyanation).
    """
    # Track functionalization steps
    functionalization_steps = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Patterns for different functional groups
                    aryl_halide_pattern = Chem.MolFromSmarts("c[Br,Cl,I]")
                    nitro_pattern = Chem.MolFromSmarts("c[N+](=[O])[O-]")
                    cyano_pattern = Chem.MolFromSmarts("cC#N")

                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants_smiles if Chem.MolFromSmiles(r)
                    ]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if product_mol and reactant_mols:
                        # Check for nitration
                        if any(
                            r.HasSubstructMatch(aryl_halide_pattern) for r in reactant_mols
                        ) and product_mol.HasSubstructMatch(nitro_pattern):
                            functionalization_steps.append(("nitration", depth))
                            print(f"Aryl halide nitration detected at depth {depth}")

                        # Check for cyanation
                        if any(
                            r.HasSubstructMatch(aryl_halide_pattern) for r in reactant_mols
                        ) and product_mol.HasSubstructMatch(cyano_pattern):
                            functionalization_steps.append(("cyanation", depth))
                            print(f"Aryl halide cyanation detected at depth {depth}")
                except Exception as e:
                    print(f"Error in aryl halide functionalization detection: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have sequential functionalization (nitration followed by cyanation)
    if len(functionalization_steps) >= 2:
        # Sort by depth (higher depth = earlier in synthesis)
        functionalization_steps.sort(key=lambda x: x[1], reverse=True)

        # Check for the sequence: nitration -> cyanation
        for i in range(len(functionalization_steps) - 1):
            if (
                functionalization_steps[i][0] == "nitration"
                and functionalization_steps[i + 1][0] == "cyanation"
            ):
                print(
                    "Sequential aryl halide functionalization detected: nitration followed by cyanation"
                )
                return True

    return False
