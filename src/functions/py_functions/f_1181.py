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
    This function detects if the synthetic route involves reduction of a nitro group to an amine.
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

            # Check if reactants contain a nitro group and product contains an amine
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and reactant_mols:
                # Check for nitro in reactants
                nitro_pattern = Chem.MolFromSmarts("[NX3+](=[OX1])[OX1-]")
                has_nitro = any(
                    mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols if mol
                )

                # Check for amine in product
                amine_pattern = Chem.MolFromSmarts("[NX3;H2]")
                has_amine = (
                    product_mol.HasSubstructMatch(amine_pattern)
                    if product_mol
                    else False
                )

                if has_nitro and has_amine:
                    print("Nitro reduction detected")
                    nitro_reduction_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_reduction_found
