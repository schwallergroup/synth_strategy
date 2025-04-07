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
    This function detects a strategy involving amide bond formation between
    an acyl chloride and an amine.
    """
    # Track if we find the pattern
    found_pattern = False

    # SMARTS patterns
    acyl_chloride = Chem.MolFromSmarts("[#6]C(=O)Cl")
    amine = Chem.MolFromSmarts("[NH2][c]")
    amide = Chem.MolFromSmarts("[#6]C(=O)[NH][c]")

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for reactants with acyl chloride and amine
            has_acyl_chloride = False
            has_amine = False

            for reactant_smiles in reactants_smiles:
                try:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(acyl_chloride):
                            has_acyl_chloride = True
                        if reactant_mol.HasSubstructMatch(amine):
                            has_amine = True
                except:
                    continue

            # Check if product has amide
            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(amide):
                    if has_acyl_chloride and has_amine:
                        print(f"Found amide formation at depth {depth}")
                        found_pattern = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_pattern
