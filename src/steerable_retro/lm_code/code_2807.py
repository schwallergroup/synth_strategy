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
    Detects if the synthetic route contains a reductive amination step
    (aldehyde + amine -> alkylated amine)
    """
    found_reductive_amination = False

    def dfs_traverse(node):
        nonlocal found_reductive_amination

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if we have an aldehyde in reactants
            aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")

            # Check if we have a secondary amine in reactants
            amine_pattern = Chem.MolFromSmarts("[N;H1]")

            # Check for tertiary amine in product (specifically with a CH2 connected)
            tertiary_amine_pattern = Chem.MolFromSmarts("[#6]-[N](-[#6])-[#6]")

            reactants = reactants_part.split(".")
            has_aldehyde = False
            has_amine = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(aldehyde_pattern):
                        has_aldehyde = True
                    if mol and mol.HasSubstructMatch(amine_pattern):
                        has_amine = True
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product_part)
                has_tertiary_amine = product_mol and product_mol.HasSubstructMatch(
                    tertiary_amine_pattern
                )
            except:
                has_tertiary_amine = False

            if has_aldehyde and has_amine and has_tertiary_amine:
                print("Found reductive amination: aldehyde + amine -> tertiary amine")
                found_reductive_amination = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_reductive_amination
