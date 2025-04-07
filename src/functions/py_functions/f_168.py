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
    Detects if the synthetic route involves nitro reduction to amine.
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
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            # Check if reactant has nitro group
            reactant_mol = Chem.MolFromSmiles(reactants_str)
            if reactant_mol:
                nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                has_nitro = reactant_mol.HasSubstructMatch(nitro_pattern)

                # Check if product has amine group where nitro was
                product_mol = Chem.MolFromSmiles(product_str)
                if product_mol:
                    amine_pattern = Chem.MolFromSmarts("[NH2]")
                    has_amine = product_mol.HasSubstructMatch(amine_pattern)

                    if (
                        has_nitro
                        and has_amine
                        and not product_mol.HasSubstructMatch(nitro_pattern)
                    ):
                        print("Detected nitro reduction to amine")
                        nitro_reduction_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_reduction_detected
