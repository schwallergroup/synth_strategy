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
    This function detects if the synthetic route employs a strategy of forming
    a heterocyclic ring (specifically pyrazole) in the early stages of synthesis.
    """
    # Track if we found a heterocycle formation
    heterocycle_formation_found = False
    heterocycle_formation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_found, heterocycle_formation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and any(reactants):
                # Check for hydrazine in reactants
                hydrazine_pattern = Chem.MolFromSmarts("[NH2][NH]")
                reactants_with_hydrazine = any(
                    mol.HasSubstructMatch(hydrazine_pattern) for mol in reactants if mol
                )

                # Check for pyrazole in product
                pyrazole_pattern = Chem.MolFromSmarts("c1nn[c]c1")  # Basic pyrazole pattern
                product_has_pyrazole = (
                    product.HasSubstructMatch(pyrazole_pattern) if product else False
                )

                # Check if this reaction forms a pyrazole
                if (
                    reactants_with_hydrazine and product_has_pyrazole and depth >= 3
                ):  # Depth >= 3 means early stage
                    print(f"Found pyrazole formation at depth {depth}")
                    heterocycle_formation_found = True
                    heterocycle_formation_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we found a heterocycle formation at early stage
    strategy_present = heterocycle_formation_found

    print(f"Heterocycle formation strategy detected: {strategy_present}")
    if strategy_present:
        print(f"Heterocycle formation occurred at depth {heterocycle_formation_depth}")

    return strategy_present
