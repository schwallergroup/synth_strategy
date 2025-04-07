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
    Detects a convergent synthesis strategy with N-alkylation as the key coupling step.
    Looks for a reaction where an amine (particularly aniline) is coupled with a brominated carbon.
    """
    found_n_alkylation = False
    fragment_count = 0

    def dfs_traverse(node):
        nonlocal found_n_alkylation, fragment_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a coupling reaction with multiple fragments
            if len(reactants) >= 2:
                fragment_count = max(fragment_count, len(reactants))

                # Check for N-alkylation pattern
                # Look for amine in one fragment and bromide in another
                has_amine = False
                has_bromide = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Check for amine/aniline pattern
                        if mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[c][NH2]")
                        ) or mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
                            has_amine = True

                        # Check for bromide pattern (typically for alkylation)
                        if mol.HasSubstructMatch(Chem.MolFromSmarts("[C][Br]")):
                            has_bromide = True

                # Check if product has the N-alkylated pattern
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[c][N][C]")):
                    if has_amine and has_bromide:
                        found_n_alkylation = True
                        print("Found N-alkylation coupling reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found a convergent N-alkylation strategy
    return found_n_alkylation and fragment_count >= 2
