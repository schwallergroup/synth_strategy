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
    This function detects thioether (C-S-C) bond formation as a fragment coupling strategy.
    """
    thioether_coupling_found = False

    def dfs_traverse(node, depth=0):
        nonlocal thioether_coupling_found

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for thioether formation
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if (
                    product_mol and len(reactant_mols) >= 2
                ):  # Need at least 2 fragments to couple
                    # Check for thioether pattern in product
                    thioether_pattern = Chem.MolFromSmarts("[c,C][S][c,C]")
                    if product_mol.HasSubstructMatch(thioether_pattern):
                        # Verify it's a new bond by checking if the pattern exists in any single reactant
                        is_new_coupling = True
                        for r_mol in reactant_mols:
                            if r_mol and r_mol.HasSubstructMatch(thioether_pattern):
                                is_new_coupling = False
                                break

                        if is_new_coupling:
                            # Check if one reactant has [S][H] or [S][c,C] and another has [Cl][c,C]
                            thiol_pattern = Chem.MolFromSmarts("[S][H]")
                            sulfide_pattern = Chem.MolFromSmarts("[S][c,C]")
                            chloro_pattern = Chem.MolFromSmarts("[Cl][c,C]")

                            has_thiol_or_sulfide = False
                            has_chloro = False

                            for r_mol in reactant_mols:
                                if r_mol:
                                    if r_mol.HasSubstructMatch(
                                        thiol_pattern
                                    ) or r_mol.HasSubstructMatch(sulfide_pattern):
                                        has_thiol_or_sulfide = True
                                    if r_mol.HasSubstructMatch(chloro_pattern):
                                        has_chloro = True

                            if has_thiol_or_sulfide and has_chloro:
                                thioether_coupling_found = True
                                print(
                                    f"Thioether fragment coupling found at depth {depth}"
                                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return thioether_coupling_found
