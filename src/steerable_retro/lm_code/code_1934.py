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
    This function detects if the synthetic route contains an ester hydrolysis
    followed by amide formation sequence.
    """
    hydrolysis_depths = []
    amide_formation_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                reactants = Chem.MolFromSmiles(reactants_smiles)
                product = Chem.MolFromSmiles(product_smiles)

                if reactants and product:
                    # Ester pattern
                    ester_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[O])")
                    # Carboxylic acid pattern
                    acid_pattern = Chem.MolFromSmarts("[#8]-[#6](=[O])")
                    # Amide pattern
                    amide_pattern = Chem.MolFromSmarts("[N;!$(N=*)][C](=[O])")
                    # Amine pattern
                    amine_pattern = Chem.MolFromSmarts("[N;!$(N=*);!$(NC=O)]")

                    # Check for ester hydrolysis
                    if (
                        reactants.HasSubstructMatch(ester_pattern)
                        and product.HasSubstructMatch(acid_pattern)
                        and not product.HasSubstructMatch(ester_pattern)
                    ):
                        hydrolysis_depths.append(depth)
                        print(f"Ester hydrolysis detected at depth {depth}: {rsmi}")

                    # Check for amide formation
                    if (
                        reactants.HasSubstructMatch(acid_pattern)
                        and reactants.HasSubstructMatch(amine_pattern)
                        and product.HasSubstructMatch(amide_pattern)
                    ):
                        amide_formation_depths.append(depth)
                        print(f"Amide formation detected at depth {depth}: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if there's at least one hydrolysis followed by amide formation
    for h_depth in hydrolysis_depths:
        for a_depth in amide_formation_depths:
            if a_depth < h_depth:  # Remember: lower depth = later in synthesis
                return True

    return False
