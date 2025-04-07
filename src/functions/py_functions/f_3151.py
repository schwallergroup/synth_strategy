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
    This function detects a nitro reduction to amine transformation in the synthesis route.
    """
    nitro_to_amine_found = False

    def dfs_traverse(node):
        nonlocal nitro_to_amine_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Create RDKit mol objects
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if all(reactant_mols) and product_mol:
                # Check for nitro group in reactants
                nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                reactants_with_nitro = any(
                    mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols
                )

                # Check for amine group in product
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                product_has_amine = product_mol.HasSubstructMatch(amine_pattern)

                # Check if nitro was reduced to amine
                if (
                    reactants_with_nitro
                    and product_has_amine
                    and not product_mol.HasSubstructMatch(nitro_pattern)
                ):
                    print("Detected nitro reduction to amine")
                    nitro_to_amine_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return nitro_to_amine_found
