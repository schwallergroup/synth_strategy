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
    This function detects if the synthesis involves multiple sequential modifications
    of an aromatic ring.
    """
    aromatic_modifications = 0

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_modifications

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Patterns for common aromatic modifications
            patterns = [
                Chem.MolFromSmarts("c[Br,Cl,I,F]"),  # Halogenation
                Chem.MolFromSmarts("c[NH2]"),  # Amination
                Chem.MolFromSmarts("c[OH]"),  # Hydroxylation
                Chem.MolFromSmarts("cO[#6]"),  # Etherification
                Chem.MolFromSmarts("c[N+](=[O])[O-]"),  # Nitration
            ]

            # Check if the reaction involves modification of aromatic ring
            prod_mol = Chem.MolFromSmiles(product)

            for reactant in reactants:
                react_mol = Chem.MolFromSmiles(reactant)

                if react_mol and prod_mol:
                    for pattern in patterns:
                        # If pattern exists in product but not in reactant, or vice versa
                        if prod_mol.HasSubstructMatch(
                            pattern
                        ) != react_mol.HasSubstructMatch(pattern):
                            aromatic_modifications += 1
                            print(f"Aromatic modification detected at depth {depth}")
                            break

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return aromatic_modifications >= 2  # At least 2 modifications
