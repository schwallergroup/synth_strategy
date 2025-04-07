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
    nitro_to_amine_found = False

    def dfs_traverse(node):
        nonlocal nitro_to_amine_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactant contains nitro group and product contains amine at same position
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and any(reactant_mol for reactant_mol in reactant_mols):
                nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")
                amine_pattern = Chem.MolFromSmarts("[#7H2]")

                for reactant_mol in reactant_mols:
                    if (
                        reactant_mol
                        and reactant_mol.HasSubstructMatch(nitro_pattern)
                        and product_mol.HasSubstructMatch(amine_pattern)
                    ):
                        print("Found nitro reduction to amine")
                        nitro_to_amine_found = True
                        break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_to_amine_found
