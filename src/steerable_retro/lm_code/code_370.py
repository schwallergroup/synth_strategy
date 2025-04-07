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
    This function detects if the synthesis route involves a sequence of
    ester reduction followed by alcohol oxidation.
    """
    ester_reduction_depths = []
    alcohol_oxidation_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester reduction
                ester_pattern = Chem.MolFromSmarts("[C](=O)[O][C]")
                alcohol_pattern = Chem.MolFromSmarts("[CH2][OH]")

                # Check for alcohol oxidation
                aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")

                # Check for ester reduction
                ester_in_reactants = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(ester_pattern)
                    for r in reactants
                    if r
                )
                alcohol_in_product = Chem.MolFromSmiles(product) is not None and Chem.MolFromSmiles(
                    product
                ).HasSubstructMatch(alcohol_pattern)

                if ester_in_reactants and alcohol_in_product:
                    print(f"Found ester reduction at depth {depth}")
                    ester_reduction_depths.append(depth)

                # Check for alcohol oxidation
                alcohol_in_reactants = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(alcohol_pattern)
                    for r in reactants
                    if r
                )
                aldehyde_in_product = Chem.MolFromSmiles(
                    product
                ) is not None and Chem.MolFromSmiles(product).HasSubstructMatch(aldehyde_pattern)

                if alcohol_in_reactants and aldehyde_in_product:
                    print(f"Found alcohol oxidation at depth {depth}")
                    alcohol_oxidation_depths.append(depth)

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have both ester reduction and alcohol oxidation in the correct sequence
    if ester_reduction_depths and alcohol_oxidation_depths:
        # In retrosynthetic direction, the alcohol oxidation should be at a lower depth
        # than the ester reduction (meaning it happens after the reduction in the forward direction)
        for er_depth in ester_reduction_depths:
            for ao_depth in alcohol_oxidation_depths:
                if (
                    ao_depth < er_depth
                ):  # In retrosynthesis, lower depth = later in forward synthesis
                    print(
                        f"Found oxidation-reduction sequence: ester reduction at depth {er_depth}, alcohol oxidation at depth {ao_depth}"
                    )
                    return True

    return False
