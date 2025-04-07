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
    This function detects the reduction of a nitro group to an amine group.
    """
    found_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal found_nitro_reduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro reduction: [N+](=O)[O-] -> [NH2]
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")

                    if reactant_mol.HasSubstructMatch(
                        nitro_pattern
                    ) and product_mol.HasSubstructMatch(amine_pattern):
                        # Ensure the nitro count decreases and amine count increases
                        nitro_count_reactant = len(reactant_mol.GetSubstructMatches(nitro_pattern))
                        nitro_count_product = len(product_mol.GetSubstructMatches(nitro_pattern))

                        amine_count_reactant = len(reactant_mol.GetSubstructMatches(amine_pattern))
                        amine_count_product = len(product_mol.GetSubstructMatches(amine_pattern))

                        if (
                            nitro_count_product < nitro_count_reactant
                            and amine_count_product > amine_count_reactant
                        ):
                            found_nitro_reduction = True
                            print("Found nitro reduction to amine")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_nitro_reduction
