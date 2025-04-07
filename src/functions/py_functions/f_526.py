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
    This function detects reductive amination for amine formation.
    """
    reductive_amination_found = False

    def dfs_traverse(node):
        nonlocal reductive_amination_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carbonyl pattern in reactants
            carbonyl_pattern = Chem.MolFromSmarts("[C](=[O])[#6]")
            # Check for amine pattern in reactants
            amine_pattern = Chem.MolFromSmarts("[N;H]")
            # Check for alkylated amine pattern in product
            alkylated_amine_pattern = Chem.MolFromSmarts("[N]([#6])[#6]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if (
                product_mol
                and any(
                    r and r.HasSubstructMatch(carbonyl_pattern) for r in reactant_mols
                )
                and any(r and r.HasSubstructMatch(amine_pattern) for r in reactant_mols)
                and product_mol.HasSubstructMatch(alkylated_amine_pattern)
            ):
                print("Reductive amination detected")
                reductive_amination_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return reductive_amination_found
