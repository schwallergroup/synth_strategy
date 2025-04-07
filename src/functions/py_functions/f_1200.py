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
    This function detects if the synthesis route employs a phenol O-alkylation with
    a morpholine-containing side chain.
    """
    morpholine_alkylation_found = False

    def dfs_traverse(node):
        nonlocal morpholine_alkylation_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phenol in reactants
            phenol_pattern = Chem.MolFromSmarts("c[OH]")

            # Check for morpholine in reactants
            morpholine_pattern = Chem.MolFromSmarts("[O]1CCN(CC)CC1")

            # Check for phenyl ether in product
            phenyl_ether_pattern = Chem.MolFromSmarts("c[O]C")

            # Check if we have at least two reactants
            if len(reactants) >= 2:
                reactant_mols = [
                    Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                ]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(phenyl_ether_pattern):
                    # Check if one reactant has phenol and another has morpholine
                    has_phenol = any(
                        mol and mol.HasSubstructMatch(phenol_pattern)
                        for mol in reactant_mols
                    )
                    has_morpholine = any(
                        mol and mol.HasSubstructMatch(morpholine_pattern)
                        for mol in reactant_mols
                    )

                    if has_phenol and has_morpholine:
                        print(
                            "Found phenol O-alkylation with morpholine-containing side chain"
                        )
                        morpholine_alkylation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return morpholine_alkylation_found
