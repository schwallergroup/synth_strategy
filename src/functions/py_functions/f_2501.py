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
    This function detects a synthetic strategy involving late-stage ether formation
    (C-O-C bond formation in the final steps).
    """
    found_late_ether_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_ether_formation

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only check final two steps (late stage)
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Count ethers in reactants and product
                prod_mol = Chem.MolFromSmiles(product)
                ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]")

                if prod_mol:
                    prod_ether_count = len(prod_mol.GetSubstructMatches(ether_pattern))

                    # Count ethers in reactants
                    reactant_ether_count = 0
                    for reactant in reactants:
                        react_mol = Chem.MolFromSmiles(reactant)
                        if react_mol:
                            reactant_ether_count += len(
                                react_mol.GetSubstructMatches(ether_pattern)
                            )

                    # If product has more ethers than reactants combined, ether formation occurred
                    if prod_ether_count > reactant_ether_count:
                        found_late_ether_formation = True
                        print(
                            f"Detected late-stage ether formation at depth {depth}: {rsmi}"
                        )
            except:
                pass

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_late_ether_formation
