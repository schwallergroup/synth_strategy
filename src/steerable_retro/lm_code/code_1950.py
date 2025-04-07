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
    This function detects a late-stage alpha-bromination of a ketone.
    """
    # Track if alpha-bromination is found
    alpha_bromination_found = False
    # Track the depth at which alpha-bromination occurs
    alpha_bromination_depth = None
    # Maximum depth in the route
    max_depth = [0]

    # Alpha-bromo ketone pattern
    alpha_bromo_ketone_pattern = Chem.MolFromSmarts("[#6]-[#6](=[O])-[#6]-[Br]")

    def dfs_traverse(node, depth=0):
        nonlocal alpha_bromination_found, alpha_bromination_depth

        # Update max depth
        max_depth[0] = max(max_depth[0], depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if reactants have alpha-bromo ketone
                reactants_have_alpha_bromo = any(
                    has_substructure(r_smi, alpha_bromo_ketone_pattern)
                    for r_smi in reactants_smiles
                )

                # Check if product has alpha-bromo ketone
                product_has_alpha_bromo = has_substructure(
                    product_smiles, alpha_bromo_ketone_pattern
                )

                # If alpha-bromo ketone is formed in this reaction
                if product_has_alpha_bromo and not reactants_have_alpha_bromo:
                    alpha_bromination_found = True
                    alpha_bromination_depth = depth
                    print(f"Alpha-bromination detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    def has_substructure(smiles, pattern):
        try:
            mol = Chem.MolFromSmiles(smiles)
            return mol and mol.HasSubstructMatch(pattern)
        except:
            return False

    # Start traversal
    dfs_traverse(route)

    # Check if alpha-bromination is in the first quarter of the synthesis (very late stage)
    if alpha_bromination_found and alpha_bromination_depth is not None:
        is_late_stage = alpha_bromination_depth <= (max_depth[0] / 4)
        print(f"Alpha-bromination is {'late-stage' if is_late_stage else 'early-stage'}")
        return is_late_stage

    return False
