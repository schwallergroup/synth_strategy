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
    This function detects if the synthetic route contains a nitro reduction step
    (converting an aromatic nitro group to an amine).
    """
    nitro_reduction_found = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if rsmi:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if reactant has nitro group and product has amine at same position
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if reactant_mol and product_mol:
                        nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
                        amine_pattern = Chem.MolFromSmarts("[NH2]")

                        if (
                            reactant_mol.HasSubstructMatch(nitro_pattern)
                            and product_mol.HasSubstructMatch(amine_pattern)
                            and not product_mol.HasSubstructMatch(nitro_pattern)
                        ):
                            print("Nitro reduction detected")
                            nitro_reduction_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return nitro_reduction_found
