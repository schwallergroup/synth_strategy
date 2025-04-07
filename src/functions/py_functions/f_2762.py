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
    This function detects if the synthetic route involves fragment coupling via ether formation.
    """
    ether_coupling_found = False

    def dfs_traverse(node):
        nonlocal ether_coupling_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Only consider reactions with multiple reactants (fragment coupling)
            if len(reactants_smiles) >= 2:
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol:
                    # Check for ether pattern in product
                    ether_pattern = Chem.MolFromSmarts("[c][OX2][#6]")

                    if product_mol.HasSubstructMatch(ether_pattern):
                        # Check if ether is not present in any of the reactants
                        ether_in_reactants = False
                        for reactant in reactants_smiles:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                ether_pattern
                            ):
                                ether_in_reactants = True
                                break

                        if not ether_in_reactants:
                            ether_coupling_found = True
                            print(
                                f"Found fragment coupling via ether formation at depth {node.get('depth', 'unknown')}"
                            )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Fragment coupling via ether formation strategy detected: {ether_coupling_found}"
    )
    return ether_coupling_found
