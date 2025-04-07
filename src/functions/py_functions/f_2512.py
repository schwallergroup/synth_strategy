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
    This function detects a strategy involving nitro group reduction.
    It looks for a reaction where a nitro group is converted to an amine.
    """
    has_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal has_nitro_reduction

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro reduction pattern
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and reactant_mols:
                # Pattern for nitro group
                nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])-[O-]")

                # Pattern for amine group
                amine_pattern = Chem.MolFromSmarts("[#6]-[N;!$(N=*);!$(N#*)]")

                # Check if any reactant has nitro group
                has_nitro_reactant = any(
                    mol and mol.HasSubstructMatch(nitro_pattern)
                    for mol in reactant_mols
                )

                # Check if product has amine group where nitro was
                has_amine_product = product_mol.HasSubstructMatch(amine_pattern)

                # Simple check: if reactant has nitro and product has amine, it might be a reduction
                if (
                    has_nitro_reactant
                    and has_amine_product
                    and not any(
                        mol and mol.HasSubstructMatch(nitro_pattern)
                        for mol in [product_mol]
                    )
                ):
                    has_nitro_reduction = True
                    print(f"Nitro reduction detected: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_nitro_reduction
