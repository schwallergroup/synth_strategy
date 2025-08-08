#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects if a nitro reduction occurs before heterocycle formation
    in the synthetic route.
    """
    # Define patterns
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    benzimidazole_pattern = Chem.MolFromSmarts("c1nc2ccccc2n1")

    # Track depths of transformations
    nitro_reduction_depth = None
    benzimidazole_formation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_depth, benzimidazole_formation_depth

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro reduction
            reactant_has_nitro = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(nitro_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            )
            product_has_nitro = (
                Chem.MolFromSmiles(product).HasSubstructMatch(nitro_pattern)
                if Chem.MolFromSmiles(product)
                else False
            )

            if reactant_has_nitro and not product_has_nitro:
                nitro_reduction_depth = depth

            # Check for benzimidazole formation
            reactants_have_benzimidazole = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(benzimidazole_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            )
            product_has_benzimidazole = (
                Chem.MolFromSmiles(product).HasSubstructMatch(benzimidazole_pattern)
                if Chem.MolFromSmiles(product)
                else False
            )

            if not reactants_have_benzimidazole and product_has_benzimidazole:
                benzimidazole_formation_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if nitro reduction occurs before heterocycle formation
    sequence_detected = (
        nitro_reduction_depth is not None
        and benzimidazole_formation_depth is not None
        and nitro_reduction_depth > benzimidazole_formation_depth
    )

    print(f"Nitro reduction detected at depth: {nitro_reduction_depth}")
    print(f"Benzimidazole formation detected at depth: {benzimidazole_formation_depth}")
    print(f"Nitro reduction before heterocycle formation: {sequence_detected}")

    return sequence_detected
