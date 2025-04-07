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
    Detects synthesis of sulfonamide-substituted piperazines, particularly
    with trifluoromethyl groups.
    """
    # Track if we found sulfonamide formation
    found_sulfonamide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_sulfonamide_formation

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for sulfonyl chloride in reactants
                sulfonyl_chloride = Chem.MolFromSmarts("[#16](=[#8])(=[#8])[Cl]")

                # Check for piperazine in reactants
                piperazine = Chem.MolFromSmarts("[#6]1[#6][#7][#6][#6][#7]1")

                # Check for sulfonamide-piperazine in product
                sulfonamide_piperazine = Chem.MolFromSmarts(
                    "[#16](=[#8])(=[#8])[#7]1[#6][#6][#7][#6][#6]1"
                )

                has_sulfonyl_chloride = False
                has_piperazine = False
                has_sulfonamide_product = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(sulfonyl_chloride):
                            has_sulfonyl_chloride = True
                        if mol.HasSubstructMatch(piperazine):
                            has_piperazine = True

                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and prod_mol.HasSubstructMatch(sulfonamide_piperazine):
                    has_sulfonamide_product = True

                if has_sulfonyl_chloride and has_piperazine and has_sulfonamide_product:
                    found_sulfonamide_formation = True
                    print(f"Found sulfonamide formation at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_sulfonamide_formation
