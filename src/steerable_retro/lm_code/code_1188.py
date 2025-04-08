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
    This function detects if the synthetic route involves nitro reduction to amine as a late-stage transformation.
    """
    late_stage_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal late_stage_nitro_reduction

        if node["type"] == "reaction" and node.get("depth", 0) == 0:  # Check if it's the final step
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Look for nitro group in reactants and amine in product
                    nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")

                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(nitro_pattern):
                                print(f"Found late-stage nitro reduction: {rsmi}")
                                late_stage_nitro_reduction = True
                except Exception as e:
                    print(f"Error in nitro reduction detection: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return late_stage_nitro_reduction
