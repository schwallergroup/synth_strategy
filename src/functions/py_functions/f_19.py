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
    This function detects a synthetic strategy where a thiazole ring is formed
    in the final step from a thioamide intermediate.
    """
    # Track if we found the pattern
    found_thiazole_formation = False
    found_thioamide_intermediate = False

    def dfs_traverse(node):
        nonlocal found_thiazole_formation, found_thioamide_intermediate

        if node["type"] == "reaction":
            # Check if this is a thiazole formation reaction (depth 0)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains thiazole
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    thiazole_pattern = Chem.MolFromSmarts("c1nc(*)sc1")
                    if product_mol.HasSubstructMatch(thiazole_pattern):
                        # Check if reactants contain thioamide
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                thioamide_pattern = Chem.MolFromSmarts("[*]C(=S)[NH2]")
                                if reactant_mol.HasSubstructMatch(thioamide_pattern):
                                    # Check if other reactant is alpha-halo carbonyl
                                    for other_reactant in reactants:
                                        if other_reactant != reactant:
                                            other_mol = Chem.MolFromSmiles(
                                                other_reactant
                                            )
                                            if other_mol:
                                                alpha_halo_carbonyl = (
                                                    Chem.MolFromSmarts(
                                                        "[*]C(=O)C[Cl,Br,I]"
                                                    )
                                                )
                                                if other_mol.HasSubstructMatch(
                                                    alpha_halo_carbonyl
                                                ):
                                                    found_thiazole_formation = True
                                                    print(
                                                        "Found thiazole formation reaction"
                                                    )

            # Check if this is a thioamide formation reaction (depth 1)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check if product contains thioamide
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    thioamide_pattern = Chem.MolFromSmarts("[*]C(=S)[NH2]")
                    if product_mol.HasSubstructMatch(thioamide_pattern):
                        found_thioamide_intermediate = True
                        print("Found thioamide intermediate")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both conditions are met
    return found_thiazole_formation and found_thioamide_intermediate
