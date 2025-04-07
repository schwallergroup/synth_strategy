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
    Detects if the synthesis route involves protection of an alcohol with a silyl group (TBS).
    """
    silyl_protection_found = False

    def dfs_traverse(node):
        nonlocal silyl_protection_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for silyl protection pattern
            # Look for alcohol in reactants and silyl ether in product
            product_mol = Chem.MolFromSmiles(product_smiles)
            silyl_ether_pattern = Chem.MolFromSmarts("[O]-[Si]([C])([C])[C]([C])([C])[C]")

            alcohol_found = False
            silyl_reagent_found = False

            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    # Check if reactant is an alcohol
                    alcohol_pattern = Chem.MolFromSmarts("[O;!H0]")
                    if reactant_mol.HasSubstructMatch(alcohol_pattern):
                        alcohol_found = True

                    # Check if reactant is a silyl chloride
                    silyl_cl_pattern = Chem.MolFromSmarts("[Cl]-[Si]")
                    if reactant_mol.HasSubstructMatch(silyl_cl_pattern):
                        silyl_reagent_found = True

            # If we found alcohol in reactants, silyl reagent, and silyl ether in product
            if (
                alcohol_found
                and silyl_reagent_found
                and product_mol
                and product_mol.HasSubstructMatch(silyl_ether_pattern)
            ):
                silyl_protection_found = True
                print("Silyl protection reaction detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return silyl_protection_found
