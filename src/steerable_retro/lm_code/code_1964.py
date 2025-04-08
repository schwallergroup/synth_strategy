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
    This function detects C1 homologation of aromatic carboxylic acid via nitrile intermediate.
    The sequence involves: ArCOOH → ArCH2OH → ArCH2X → ArCH2CN → ArCH2COOH → ArCH2COOMe
    """
    # Initialize tracking variables
    has_aromatic_acid = False
    has_benzylic_alcohol = False
    has_benzylic_halide = False
    has_phenylacetonitrile = False
    has_phenylacetic_acid = False
    has_phenylacetic_ester = False

    # SMARTS patterns
    aromatic_acid_pattern = Chem.MolFromSmarts("[c][C](=O)[OH]")
    benzylic_alcohol_pattern = Chem.MolFromSmarts("[c][CH2][OH]")
    benzylic_halide_pattern = Chem.MolFromSmarts("[c][CH2][Br,Cl,I,F]")
    phenylacetonitrile_pattern = Chem.MolFromSmarts("[c][CH2][C]#[N]")
    phenylacetic_acid_pattern = Chem.MolFromSmarts("[c][CH2][C](=O)[OH]")
    phenylacetic_ester_pattern = Chem.MolFromSmarts("[c][CH2][C](=O)[O][C]")

    def dfs_traverse(node):
        nonlocal has_aromatic_acid, has_benzylic_alcohol, has_benzylic_halide
        nonlocal has_phenylacetonitrile, has_phenylacetic_acid, has_phenylacetic_ester

        if node["type"] == "mol":
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    if mol.HasSubstructMatch(aromatic_acid_pattern):
                        has_aromatic_acid = True
                    if mol.HasSubstructMatch(benzylic_alcohol_pattern):
                        has_benzylic_alcohol = True
                    if mol.HasSubstructMatch(benzylic_halide_pattern):
                        has_benzylic_halide = True
                    if mol.HasSubstructMatch(phenylacetonitrile_pattern):
                        has_phenylacetonitrile = True
                    if mol.HasSubstructMatch(phenylacetic_acid_pattern):
                        has_phenylacetic_acid = True
                    if mol.HasSubstructMatch(phenylacetic_ester_pattern):
                        has_phenylacetic_ester = True
            except:
                print(f"Error processing molecule: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Process reaction nodes if needed
            pass

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have the complete homologation sequence
    sequence_present = (
        has_aromatic_acid
        and has_benzylic_alcohol
        and has_benzylic_halide
        and has_phenylacetonitrile
        and has_phenylacetic_acid
    )

    print(f"C1 homologation via nitrile detection results:")
    print(f"  Aromatic acid: {has_aromatic_acid}")
    print(f"  Benzylic alcohol: {has_benzylic_alcohol}")
    print(f"  Benzylic halide: {has_benzylic_halide}")
    print(f"  Phenylacetonitrile: {has_phenylacetonitrile}")
    print(f"  Phenylacetic acid: {has_phenylacetic_acid}")
    print(f"  Phenylacetic ester: {has_phenylacetic_ester}")
    print(f"  Complete sequence present: {sequence_present}")

    return sequence_present
