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
    Detects the specific nitrogen transformation sequence:
    nitro → indole NH → tertiary amine → nitrile
    """
    nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])-[O-]")
    indole_nh_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")
    tertiary_amine_pattern = Chem.MolFromSmarts("[#6]-[N](-[#6])-[#6]")
    nitrile_pattern = Chem.MolFromSmarts("[#6]#[N]")

    # Track the sequence of nitrogen transformations
    transformation_sequence = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product)
                if not product_mol:
                    return

                # Check for nitrogen functional groups in product
                if product_mol.HasSubstructMatch(nitro_pattern):
                    transformation_sequence.append(("nitro", depth))
                if product_mol.HasSubstructMatch(indole_nh_pattern):
                    transformation_sequence.append(("indole_nh", depth))
                if product_mol.HasSubstructMatch(tertiary_amine_pattern):
                    transformation_sequence.append(("tertiary_amine", depth))
                if product_mol.HasSubstructMatch(nitrile_pattern):
                    transformation_sequence.append(("nitrile", depth))

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Sort by depth (descending) to get chronological order in synthesis direction
    transformation_sequence.sort(key=lambda x: x[1], reverse=True)

    # Extract just the functional group names
    sequence = [item[0] for item in transformation_sequence]
    print(f"Nitrogen transformation sequence: {sequence}")

    # Check if our target sequence is a subsequence of the detected sequence
    target_sequence = ["nitro", "indole_nh", "tertiary_amine", "nitrile"]

    # Check if target_sequence appears in sequence in order (not necessarily consecutive)
    i, j = 0, 0
    while i < len(sequence) and j < len(target_sequence):
        if sequence[i] == target_sequence[j]:
            j += 1
        i += 1

    return j == len(target_sequence)
