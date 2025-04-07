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
    This function detects if the synthetic route involves the formation of ether linkages,
    particularly aryl ethers (Ar-O-C).
    """
    ether_formation_found = False

    def dfs_traverse(node):
        nonlocal ether_formation_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains aryl ether
            if product:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Check for aryl ether pattern
                    aryl_ether_pattern = Chem.MolFromSmarts("[c]-[O]-[C]")
                    has_aryl_ether_product = product_mol.HasSubstructMatch(
                        aryl_ether_pattern
                    )

                    # Check if reactants contain the same pattern
                    has_aryl_ether_reactant = False

                    for reactant in reactants:
                        if reactant:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                aryl_ether_pattern
                            ):
                                has_aryl_ether_reactant = True
                                break

                    # If product has aryl ether that's not in reactants, it's formed in this reaction
                    if has_aryl_ether_product and not has_aryl_ether_reactant:
                        print("Aryl ether linkage formation detected")
                        ether_formation_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return ether_formation_found
