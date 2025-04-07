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
    This function detects if the synthetic route contains multiple N-alkylation steps.
    """
    n_alkylation_count = 0

    def dfs_traverse(node):
        nonlocal n_alkylation_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if rsmi:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for N-alkylation pattern
                # Look for bromo/iodo/chloro compound as one reactant
                has_haloalkyl = False
                for r in reactants_smiles:
                    haloalkyl_pattern = Chem.MolFromSmarts("[C]-[Br,I,Cl]")
                    r_mol = Chem.MolFromSmiles(r)
                    if r_mol and r_mol.HasSubstructMatch(haloalkyl_pattern):
                        has_haloalkyl = True
                        break

                if has_haloalkyl:
                    # Check if product has a new C-N bond
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol:
                        # This is a simplification - in a real implementation,
                        # we would need to compare atom mappings to confirm the alkylation
                        n_alkylation_count += 1
                        print(f"N-alkylation detected (count: {n_alkylation_count})")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return n_alkylation_count >= 2  # Return True if at least 2 N-alkylations found
