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
    Detects reductive amination strategy (aldehyde/ketone + amine → secondary/tertiary amine).
    """
    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction":
            # Check if this is a reaction node
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for reductive amination patterns
                aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
                amine_pattern = Chem.MolFromSmarts("[NH1,NH2]")
                sec_amine_pattern = Chem.MolFromSmarts("[NH1][CH2]")

                # Check reactants for aldehyde and amine
                has_aldehyde = False
                has_amine = False

                for r in reactants:
                    try:
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            if mol.HasSubstructMatch(aldehyde_pattern):
                                has_aldehyde = True
                            if mol.HasSubstructMatch(amine_pattern):
                                has_amine = True
                    except:
                        continue

                # Check product for secondary amine
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    has_sec_amine = prod_mol and prod_mol.HasSubstructMatch(sec_amine_pattern)
                except:
                    has_sec_amine = False

                # If we have aldehyde, amine, and secondary amine in product, it's likely reductive amination
                if has_aldehyde and has_amine and has_sec_amine:
                    found_pattern = True
                    print(f"Found reductive amination at depth {depth}")

        # Traverse children
        if "children" in node:
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Reductive amination strategy: {found_pattern}")
    return found_pattern
