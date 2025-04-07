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
    This function detects the reduction of a nitro group to an amine
    in the synthetic route.
    """
    has_nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal has_nitro_reduction

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for nitro reduction
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol:
                    # Nitro pattern and amine pattern
                    nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")

                    # Check if any reactant has nitro group
                    reactant_has_nitro = any(
                        mol and mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols
                    )

                    # Check if product has amine group
                    product_has_amine = product_mol.HasSubstructMatch(amine_pattern)

                    # If reactant has nitro and product has amine, nitro reduction occurred
                    if reactant_has_nitro and product_has_amine:
                        print(f"Found nitro to amine reduction at depth {depth}")
                        has_nitro_reduction = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_nitro_reduction
