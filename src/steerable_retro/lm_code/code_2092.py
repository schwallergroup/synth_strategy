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
    This function detects if the synthesis follows a sequential functionalization pattern
    after heterocycle formation, including N-alkylation, reduction, and ether formation.
    """
    # Track key transformations
    n_alkylation_detected = False
    reduction_detected = False
    ether_formation_detected = False

    def dfs_traverse(node):
        nonlocal n_alkylation_detected, reduction_detected, ether_formation_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product_smiles)

            # Check for N-alkylation (N-H to N-C)
            if len(reactants_smiles) == 2:
                for r_smiles in reactants_smiles:
                    r_mol = Chem.MolFromSmiles(r_smiles)
                    if r_mol and product_mol:
                        # Check for benzimidazole N-H
                        nh_pattern = Chem.MolFromSmarts("[nH]1cnc2ccccc12")
                        # Check for N-alkylated product
                        n_alkyl_pattern = Chem.MolFromSmarts("[n]1([C])cnc2ccccc12")

                        if r_mol.HasSubstructMatch(nh_pattern) and product_mol.HasSubstructMatch(
                            n_alkyl_pattern
                        ):
                            n_alkylation_detected = True
                            print("N-alkylation step detected")

            # Check for ester reduction (C(=O)OC to CH2OH)
            ester_pattern = Chem.MolFromSmarts("C(=O)OC")
            alcohol_pattern = Chem.MolFromSmarts("CO")

            for r_smiles in reactants_smiles:
                r_mol = Chem.MolFromSmiles(r_smiles)
                if r_mol and product_mol:
                    if r_mol.HasSubstructMatch(ester_pattern) and product_mol.HasSubstructMatch(
                        alcohol_pattern
                    ):
                        # Additional check to confirm reduction
                        if not product_mol.HasSubstructMatch(ester_pattern):
                            reduction_detected = True
                            print("Ester reduction step detected")

            # Check for ether formation
            for r_smiles in reactants_smiles:
                r_mol = Chem.MolFromSmiles(r_smiles)
                if r_mol and product_mol:
                    # Check for alcohol in reactant
                    alcohol_pattern = Chem.MolFromSmarts("CO")
                    # Check for ether in product
                    ether_pattern = Chem.MolFromSmarts("COC")

                    if r_mol.HasSubstructMatch(alcohol_pattern) and product_mol.HasSubstructMatch(
                        ether_pattern
                    ):
                        ether_formation_detected = True
                        print("Ether formation step detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have the sequential functionalization pattern
    has_sequential_functionalization = n_alkylation_detected and (
        reduction_detected or ether_formation_detected
    )

    print(f"Sequential functionalization strategy detected: {has_sequential_functionalization}")
    return has_sequential_functionalization
