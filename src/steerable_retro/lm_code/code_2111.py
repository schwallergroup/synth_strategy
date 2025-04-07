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
    This function detects if the synthesis route involves a late-stage N-oxidation
    (N-oxidation in the final or penultimate step).
    """
    n_oxidation_detected = False
    depth_of_n_oxidation = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal n_oxidation_detected, depth_of_n_oxidation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for N-oxidation: presence of N+-O- in product but not in reactants
            product_mol = Chem.MolFromSmiles(product_smiles)
            n_oxide_pattern = Chem.MolFromSmarts("[#7+]-[O-]")

            if product_mol and product_mol.HasSubstructMatch(n_oxide_pattern):
                # Check if N-oxide was not present in reactants
                n_oxide_in_reactants = False
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(n_oxide_pattern):
                        n_oxide_in_reactants = True
                        break

                if not n_oxide_in_reactants:
                    n_oxidation_detected = True
                    depth_of_n_oxidation = depth
                    print(f"N-oxidation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if N-oxidation occurred in the first two steps (late stage)
    # Remember: depth 0 or 1 means late in the synthesis
    is_late_stage = depth_of_n_oxidation <= 1

    if n_oxidation_detected and is_late_stage:
        print("Late-stage N-oxidation strategy detected")
        return True
    return False
