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
    This function detects nitro reduction to amine in the synthetic route.
    """
    nitro_reduction_detected = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Try to detect nitro reduction pattern
            try:
                p_mol = Chem.MolFromSmiles(product)
                amine_pattern = Chem.MolFromSmarts("[N;H2]-[c,C]")

                for r in reactants:
                    r_mol = Chem.MolFromSmiles(r)
                    if r_mol:
                        nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
                        if r_mol.HasSubstructMatch(
                            nitro_pattern
                        ) and p_mol.HasSubstructMatch(amine_pattern):
                            print("Nitro reduction detected in reaction:", rsmi)
                            nitro_reduction_detected = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_reduction_detected
