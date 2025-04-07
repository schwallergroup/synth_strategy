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
    Detects if the synthesis route includes an O-methylation step.
    """
    o_methylation_found = False

    def dfs_traverse(node):
        nonlocal o_methylation_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check if reactants contain an alcohol
            reactants = reactants_smiles.split(".")
            has_alcohol = any(
                Chem.MolFromSmiles(r)
                and Chem.MolFromSmiles(r).HasSubstructMatch(Chem.MolFromSmarts("[#8;H1]"))
                for r in reactants
                if r
            )

            # Check if product contains a methyl ether
            methyl_ether_pattern = Chem.MolFromSmarts("[#8][#6;H3]")
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None
            has_methyl_ether = product_mol and product_mol.HasSubstructMatch(methyl_ether_pattern)

            if has_alcohol and has_methyl_ether:
                print("Found O-methylation step")
                o_methylation_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return o_methylation_found
