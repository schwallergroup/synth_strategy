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
    Detects if the synthesis involves reduction of a nitro group to an amine.
    """
    has_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal has_nitro_reduction

        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Create molecules
            product_mol = Chem.MolFromSmiles(product)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

            if product_mol and all(reactant_mols):
                # Check for nitro reduction pattern
                nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")
                amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")

                has_nitro = False
                for r_mol in reactant_mols:
                    if r_mol.HasSubstructMatch(nitro_pattern):
                        has_nitro = True
                        break

                has_amine = product_mol.HasSubstructMatch(amine_pattern)

                # If we have a nitro group as reactant and an amine as product,
                # it's likely a nitro reduction
                if has_nitro and has_amine:
                    has_nitro_reduction = True
                    print(f"Nitro reduction detected: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_nitro_reduction
