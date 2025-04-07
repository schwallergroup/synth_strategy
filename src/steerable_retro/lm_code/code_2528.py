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
    This function detects a strategy involving late-stage amide coupling with a piperazine fragment.
    """
    # Track if we found the required transformations
    amide_coupling_with_piperazine = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_with_piperazine

        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions (low depth)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Convert to RDKit molecules
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product = Chem.MolFromSmiles(product_smiles)

                    if product and all(r for r in reactants):
                        # Check for amide formation
                        amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
                        carboxylic_pattern = Chem.MolFromSmarts("[C](=[O])[O;H]")
                        piperazine_pattern = Chem.MolFromSmarts("[N]1[C][C][N][C][C]1")

                        # Check if product has amide bond and piperazine
                        if (
                            product.HasSubstructMatch(amide_pattern)
                            and product.HasSubstructMatch(piperazine_pattern)
                            and any(r.HasSubstructMatch(carboxylic_pattern) for r in reactants)
                        ):
                            amide_coupling_with_piperazine = True
                            print("Detected late-stage amide coupling with piperazine")
                except:
                    print("Error processing reaction SMILES")

        # Process children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return amide_coupling_with_piperazine
