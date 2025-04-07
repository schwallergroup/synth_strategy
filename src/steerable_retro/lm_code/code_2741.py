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
    This function detects a synthetic strategy involving reduction of a nitro group to an amine.
    """
    nitro_reduction_detected = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if all(reactants_mols) and product_mol:
                    # Check for nitro group in reactants
                    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                    reactants_with_nitro = any(
                        [mol.HasSubstructMatch(nitro_pattern) for mol in reactants_mols]
                    )

                    # Check for amine in product
                    amine_pattern = Chem.MolFromSmarts("[NH2]")
                    product_has_amine = product_mol.HasSubstructMatch(amine_pattern)

                    if reactants_with_nitro and product_has_amine:
                        print("Nitro to amine reduction detected")
                        nitro_reduction_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return nitro_reduction_detected
