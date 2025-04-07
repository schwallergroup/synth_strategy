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
    Detects if the synthesis route includes a dibenzylation step.
    """
    dibenzylation_found = False

    def dfs_traverse(node):
        nonlocal dibenzylation_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check if reactants contain a primary amine
            reactants = reactants_smiles.split(".")
            primary_amine_pattern = Chem.MolFromSmarts("[#7;H2]")
            has_primary_amine = any(
                Chem.MolFromSmiles(r)
                and Chem.MolFromSmiles(r).HasSubstructMatch(primary_amine_pattern)
                for r in reactants
                if r
            )

            # Check if product contains a dibenzylated amine
            dibenzyl_pattern = Chem.MolFromSmarts("[#7]([#6][c])([#6][c])")
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None
            has_dibenzyl = product_mol and product_mol.HasSubstructMatch(
                dibenzyl_pattern
            )

            if has_primary_amine and has_dibenzyl:
                print("Found dibenzylation step")
                dibenzylation_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return dibenzylation_found
