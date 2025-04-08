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
    Detects if the synthesis route involves a late-stage coupling between
    a heterocyclic system and an aliphatic ring system.
    """
    late_stage_coupling = False
    max_depth = 2  # Consider reactions with depth <= 2 as late-stage

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_coupling

        if node["type"] == "reaction" and depth <= max_depth:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is a coupling reaction between heterocycle and aliphatic
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and len(reactant_mols) >= 2 and all(r for r in reactant_mols):
                # Check for heterocycle in reactants
                heterocycle_pattern = Chem.MolFromSmarts("[r5,r6]1:[a]:[a]:[a]:[a]:1")
                # Check for aliphatic ring in reactants
                aliphatic_pattern = Chem.MolFromSmarts(
                    "[CR1,CR2,CR3]1[CR1,CR2,CR3][CR1,CR2,CR3][CR1,CR2,CR3][CR1,CR2,CR3][CR1,CR2,CR3]1"
                )

                has_heterocycle = any(
                    mol.HasSubstructMatch(heterocycle_pattern) for mol in reactant_mols
                )
                has_aliphatic = any(
                    mol.HasSubstructMatch(aliphatic_pattern) for mol in reactant_mols
                )

                # Check if product has both heterocycle and aliphatic connected
                if (
                    has_heterocycle
                    and has_aliphatic
                    and product_mol.HasSubstructMatch(heterocycle_pattern)
                    and product_mol.HasSubstructMatch(aliphatic_pattern)
                ):
                    print(f"Found late-stage heterocycle-aliphatic coupling at depth {depth}")
                    late_stage_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_coupling
