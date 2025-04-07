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
    Detects a specific sequence of functional group interconversions:
    acid → ester → alcohol → bromide → amine
    """
    # Track the sequence of transformations
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Define patterns
            acid_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#8;H1]")
            ester_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#8]-[#6]")
            alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")
            bromide_pattern = Chem.MolFromSmarts("[#6]-[#35]")
            amine_pattern = Chem.MolFromSmarts("[#6]-[#7]")

            # Check for acid to ester
            if any(
                r.HasSubstructMatch(acid_pattern) for r in reactants if r is not None
            ) and product.HasSubstructMatch(ester_pattern):
                transformations.append(("acid_to_ester", depth))
                print(f"Acid to ester transformation at depth {depth}")

            # Check for ester to alcohol
            if any(
                r.HasSubstructMatch(ester_pattern) for r in reactants if r is not None
            ) and product.HasSubstructMatch(alcohol_pattern):
                transformations.append(("ester_to_alcohol", depth))
                print(f"Ester to alcohol transformation at depth {depth}")

            # Check for alcohol to bromide
            if any(
                r.HasSubstructMatch(alcohol_pattern) for r in reactants if r is not None
            ) and product.HasSubstructMatch(bromide_pattern):
                transformations.append(("alcohol_to_bromide", depth))
                print(f"Alcohol to bromide transformation at depth {depth}")

            # Check for bromide to amine (N-alkylation)
            if any(
                r.HasSubstructMatch(bromide_pattern) for r in reactants if r is not None
            ) and product.HasSubstructMatch(amine_pattern):
                transformations.append(("bromide_to_amine", depth))
                print(f"Bromide to amine transformation at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the sequence is present in the correct order
    # Note: In retrosynthetic direction, higher depths are earlier in the synthesis
    # So we need to check if the sequence appears in reverse order by depth

    # Extract transformation types in order of increasing depth (earlier in synthesis)
    ordered_transformations = sorted(transformations, key=lambda x: x[1], reverse=True)
    transformation_types = [t[0] for t in ordered_transformations]

    # Check for the sequence
    sequence = ["acid_to_ester", "ester_to_alcohol", "alcohol_to_bromide", "bromide_to_amine"]

    # Check if sequence appears as a subsequence
    def is_subsequence(subseq, seq):
        it = iter(seq)
        return all(c in it for c in subseq)

    result = is_subsequence(sequence, transformation_types)

    if result:
        print("Functional group interconversion sequence detected")
    else:
        print("Functional group interconversion sequence not detected")
        print(f"Found transformations: {transformation_types}")

    return result
