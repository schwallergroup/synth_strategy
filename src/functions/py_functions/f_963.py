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
    This function detects sequences of functional group interconversions,
    particularly focusing on transformations between alkenes, aldehydes,
    alcohols, esters, and carboxylic acids.
    """
    # Track functional group changes through the synthesis
    fg_sequence = []
    has_complex_sequence = False

    # Define SMARTS patterns for functional groups
    alkene_pattern = Chem.MolFromSmarts("[#6]=[#6]")
    aldehyde_pattern = Chem.MolFromSmarts("[#6](=[#8])-[#1]")
    alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#1]")
    ester_pattern = Chem.MolFromSmarts("[#6](=[#8])-[#8]-[#6]")
    carboxylic_acid_pattern = Chem.MolFromSmarts("[#6](=[#8])-[#8]-[#1]")

    def identify_functional_group(mol):
        """Helper function to identify the primary functional group in a molecule"""
        if not mol:
            return "unknown"

        if mol.HasSubstructMatch(carboxylic_acid_pattern):
            return "carboxylic_acid"
        if mol.HasSubstructMatch(ester_pattern):
            return "ester"
        if mol.HasSubstructMatch(aldehyde_pattern):
            return "aldehyde"
        if mol.HasSubstructMatch(alcohol_pattern):
            return "alcohol"
        if mol.HasSubstructMatch(alkene_pattern):
            return "alkene"

        return "other"

    def dfs_traverse(node, depth=0):
        nonlocal fg_sequence

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]
            product_mol = Chem.MolFromSmiles(product)

            fg = identify_functional_group(product_mol)
            fg_sequence.append((depth, fg))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Sort by depth (ascending)
    fg_sequence.sort(key=lambda x: x[0])

    # Extract just the functional groups in sequence
    fg_types = [fg for _, fg in fg_sequence]

    # Check for complex sequences (3+ different functional groups)
    unique_fgs = set(fg_types)
    if len(unique_fgs) >= 3:
        has_complex_sequence = True
        print(
            f"Complex functional group interconversion sequence detected: {' → '.join(fg_types)}"
        )

    return has_complex_sequence
