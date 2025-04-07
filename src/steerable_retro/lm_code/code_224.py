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
    This function detects a borylation of aryl bromide followed by Suzuki coupling strategy.
    It looks for conversion of aryl bromide to boronic acid/ester and subsequent biaryl formation.
    """
    # Track if we've found borylation and coupling steps
    borylation_found = False
    suzuki_coupling_found = False

    def dfs_traverse(node):
        nonlocal borylation_found, suzuki_coupling_found

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for borylation (aryl bromide to boronic acid/ester)
            if not borylation_found:
                aryl_bromide_pattern = Chem.MolFromSmarts("[c]-[Br]")
                boronic_pattern = Chem.MolFromSmarts("[c]-[B]")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if (
                    product_mol
                    and any(r and r.HasSubstructMatch(aryl_bromide_pattern) for r in reactant_mols)
                    and product_mol.HasSubstructMatch(boronic_pattern)
                ):
                    borylation_found = True
                    print("Found borylation step")

            # Check for Suzuki coupling (boronic acid + aryl halide/triflate → biaryl)
            if not suzuki_coupling_found:
                boronic_pattern = Chem.MolFromSmarts("[c]-[B]")
                biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if (
                    product_mol
                    and any(r and r.HasSubstructMatch(boronic_pattern) for r in reactant_mols)
                    and product_mol.HasSubstructMatch(biaryl_pattern)
                ):
                    suzuki_coupling_found = True
                    print("Found Suzuki coupling step")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both borylation and Suzuki coupling were found
    return borylation_found and suzuki_coupling_found
