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
    Detects a strategy involving multiple heteroatom linkages (ether and thioether)
    between aromatic rings with carboxylic acid protection-deprotection.
    """
    # Track key features
    has_thioether_formation = False
    has_diaryl_ether_formation = False
    has_carboxylic_protection = False
    has_carboxylic_deprotection = False

    def dfs_traverse(node):
        nonlocal has_thioether_formation, has_diaryl_ether_formation
        nonlocal has_carboxylic_protection, has_carboxylic_deprotection

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and len(reactants) >= 1:
                # Check for thioether formation (C-S-C between aromatics)
                thioether_pattern = Chem.MolFromSmarts("[c][CH2][S][CH2][c]")
                if product.HasSubstructMatch(thioether_pattern):
                    # Check if this is a new formation by ensuring reactants don't have it
                    if not any(
                        r.HasSubstructMatch(thioether_pattern) for r in reactants if r
                    ):
                        print("Detected thioether formation")
                        has_thioether_formation = True

                # Check for diaryl ether formation (C-O-C between aromatics)
                diaryl_ether_pattern = Chem.MolFromSmarts("[c][O][c]")
                if product.HasSubstructMatch(diaryl_ether_pattern):
                    # Check if this is a new formation by ensuring reactants don't have it
                    if not any(
                        r.HasSubstructMatch(diaryl_ether_pattern)
                        for r in reactants
                        if r
                    ):
                        print("Detected diaryl ether formation")
                        has_diaryl_ether_formation = True

                # Check for carboxylic acid protection (acid to ester)
                acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
                ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")

                if any(
                    r.HasSubstructMatch(acid_pattern) for r in reactants if r
                ) and product.HasSubstructMatch(ester_pattern):
                    print("Detected carboxylic acid protection")
                    has_carboxylic_protection = True

                # Check for carboxylic acid deprotection (ester to acid)
                if any(
                    r.HasSubstructMatch(ester_pattern) for r in reactants if r
                ) and product.HasSubstructMatch(acid_pattern):
                    print("Detected carboxylic acid deprotection")
                    has_carboxylic_deprotection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have both types of heteroatom linkages and acid protection/deprotection
    return (
        has_thioether_formation
        and has_diaryl_ether_formation
        and (has_carboxylic_protection or has_carboxylic_deprotection)
    )
