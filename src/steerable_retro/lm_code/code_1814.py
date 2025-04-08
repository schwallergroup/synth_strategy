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
    Detects formation of ether linkage between aromatic systems.
    Looks for C-O-C bond formation where both carbons are part of aromatic systems.
    """
    ether_formation_detected = False

    def dfs_traverse(node):
        nonlocal ether_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alcohol in reactants
                alcohol_pattern = Chem.MolFromSmarts("[c][C][OH]")
                # Check for phenol in reactants
                phenol_pattern = Chem.MolFromSmarts("[c][OH]")
                # Check for ether in product
                ether_pattern = Chem.MolFromSmarts("[c][C][O][c]")

                has_alcohol = False
                has_phenol = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(alcohol_pattern):
                            has_alcohol = True
                        if mol.HasSubstructMatch(phenol_pattern):
                            has_phenol = True

                prod_mol = Chem.MolFromSmiles(product)
                if (
                    prod_mol
                    and has_alcohol
                    and has_phenol
                    and prod_mol.HasSubstructMatch(ether_pattern)
                ):
                    ether_formation_detected = True
                    print(f"Detected ether formation between aromatic systems: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return ether_formation_detected
