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
    Detects if the synthetic route involves a nitro reduction to amine.
    """
    nitro_reduction_found = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant contains a nitro group
            nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            for reactant in reactants:
                try:
                    r_mol = Chem.MolFromSmiles(reactant)
                    if r_mol and r_mol.HasSubstructMatch(nitro_pattern):
                        p_mol = Chem.MolFromSmiles(product)
                        # Check if product has an amine where nitro was
                        if p_mol and p_mol.HasSubstructMatch(amine_pattern):
                            print(f"Nitro reduction detected: {rsmi}")
                            nitro_reduction_found = True
                except:
                    continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return nitro_reduction_found
