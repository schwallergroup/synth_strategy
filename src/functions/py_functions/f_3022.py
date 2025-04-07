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
    Detects if the synthesis route involves an amide formation step.
    """
    amide_formed = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_formed

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol:
                # Check for amide pattern in product
                amide_pattern = Chem.MolFromSmarts("[#7]-[#6](=[#8])")

                if product_mol.HasSubstructMatch(amide_pattern):
                    # Check if amide was not present in reactants
                    reactants_have_amide = any(
                        r and r.HasSubstructMatch(amide_pattern)
                        for r in reactant_mols
                        if r
                    )

                    if not reactants_have_amide:
                        amide_formed = True
                        print(f"Amide formation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return amide_formed
