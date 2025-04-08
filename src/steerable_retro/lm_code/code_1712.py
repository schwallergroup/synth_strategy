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
    Detects a strategy involving the reduction of a nitro group to an amine
    early in the synthesis.
    """
    found_nitro_reduction = False

    # SMARTS patterns
    nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])-[O-]")
    amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")

    def dfs_traverse(node, depth=0):
        nonlocal found_nitro_reduction

        if node["type"] == "reaction" and depth >= 3:  # Early in synthesis (high depth)
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Create RDKit molecules
            product_mol = Chem.MolFromSmiles(product)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

            # Check for nitro reduction
            if (
                any(mol and mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols)
                and product_mol
                and product_mol.HasSubstructMatch(amine_pattern)
            ):
                found_nitro_reduction = True
                print(f"Found nitro reduction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitro reduction strategy present: {found_nitro_reduction}")
    return found_nitro_reduction
