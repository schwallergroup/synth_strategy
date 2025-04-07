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
    This function detects if the synthesis follows a linear strategy
    (each step has only one non-commercial reactant).
    """
    is_linear = True

    # First, build a dictionary mapping SMILES to their nodes for quick lookup
    smiles_to_node = {}

    def build_smiles_map(node):
        if node["type"] == "mol" and "smiles" in node:
            # Canonicalize SMILES for consistent matching
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                canonical_smiles = Chem.MolToSmiles(mol)
                smiles_to_node[canonical_smiles] = node

        # Process children recursively
        for child in node.get("children", []):
            build_smiles_map(child)

    # Build the SMILES to node mapping
    build_smiles_map(route)

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if not is_linear:  # Early return if we already found it's not linear
            return

        if node["type"] == "reaction":
            try:
                if "metadata" in node and "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    reactants_part = rsmi.split(">")[0]
                    reactants = reactants_part.split(".")

                    # Count non-commercial reactants
                    non_commercial_reactants = 0

                    for r in reactants:
                        if not r:  # Skip empty strings
                            continue

                        # Canonicalize reactant SMILES
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            canonical_smiles = Chem.MolToSmiles(mol)

                            # Find the corresponding node
                            reactant_node = smiles_to_node.get(canonical_smiles)

                            # Check if it's a non-commercial reactant
                            if reactant_node and not reactant_node.get(
                                "in_stock", False
                            ):
                                non_commercial_reactants += 1

                    if non_commercial_reactants > 1:
                        print(
                            f"Found convergent step with {non_commercial_reactants} non-commercial reactants at depth {depth}"
                        )
                        is_linear = False
            except Exception as e:
                print(f"Error processing reaction node at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return is_linear
