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
    This function detects a synthetic strategy where borylation is used
    to prepare an aromatic ring for subsequent cross-coupling.
    """
    has_borylation = False

    def dfs_traverse(node):
        nonlocal has_borylation

        if node["type"] == "reaction":
            # Extract reaction SMILES
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            # Split into reactants and product
            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants = parts[0].split(".")
            product = parts[2]

            # Check for borylation patterns
            try:
                # Look for boronic acid in product but not in reactants
                boronic_acid_pattern = Chem.MolFromSmarts("[c][B]([OH])[OH]")

                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                product_has_boronic_acid = (
                    product_mol and product_mol.HasSubstructMatch(boronic_acid_pattern)
                )
                reactants_have_boronic_acid = any(
                    mol and mol.HasSubstructMatch(boronic_acid_pattern)
                    for mol in reactant_mols
                )

                if product_has_boronic_acid and not reactants_have_boronic_acid:
                    depth = node["metadata"].get("depth", -1)
                    print(f"Detected borylation preparation at depth {depth}")
                    has_borylation = True
            except:
                pass

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return has_borylation
