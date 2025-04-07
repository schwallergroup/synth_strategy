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
    This function detects a synthetic strategy where a linear molecule is transformed
    into a bicyclic system in the final steps.
    """
    bicyclic_formed = False
    bicyclic_depth = float("inf")

    def is_bicyclic(mol):
        """Check if a molecule has a true bicyclic system (rings sharing atoms)"""
        if mol is None:
            return False

        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()

        # Need at least 2 rings to be bicyclic
        if len(rings) < 2:
            return False

        # Check if any pair of rings shares at least one atom
        for i in range(len(rings)):
            for j in range(i + 1, len(rings)):
                if set(rings[i]).intersection(set(rings[j])):
                    return True

        return False

    def is_linear_or_monocyclic(mol):
        """Check if a molecule is linear (no rings) or has only one ring"""
        if mol is None:
            return False

        ring_info = mol.GetRingInfo()
        return ring_info.NumRings() <= 1

    def dfs_traverse(node, depth=0):
        nonlocal bicyclic_formed, bicyclic_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if all(reactant_mols) and product_mol:
                    # Check if product is bicyclic and all reactants are linear or monocyclic
                    if is_bicyclic(product_mol) and all(
                        is_linear_or_monocyclic(r) for r in reactant_mols
                    ):
                        print(f"Bicyclic system formation detected at depth {depth}")
                        bicyclic_formed = True
                        bicyclic_depth = min(bicyclic_depth, depth)
            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if bicyclic system was formed in the late stage (depth <= 3)
    if bicyclic_formed and bicyclic_depth <= 3:
        print(f"Late-stage bicyclic formation detected at depth {bicyclic_depth}")
        return True
    return False
