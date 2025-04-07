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
    This function detects if the synthetic route involves formation of C-N or C-O bonds.
    """
    heteroatom_bond_formation_found = False

    def dfs_traverse(node):
        nonlocal heteroatom_bond_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check if we can detect C-N or C-O bond formation
            if product:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Look for C-N and C-O bonds in product
                    c_n_pattern = Chem.MolFromSmarts("[#6]-[#7]")
                    c_o_pattern = Chem.MolFromSmarts("[#6]-[#8]")

                    c_n_matches_product = product_mol.GetSubstructMatches(c_n_pattern)
                    c_o_matches_product = product_mol.GetSubstructMatches(c_o_pattern)

                    # Check if these bonds exist in reactants
                    reactants = reactants_part.split(".")
                    c_n_matches_reactants = []
                    c_o_matches_reactants = []

                    for reactant in reactants:
                        if reactant:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                c_n_matches_reactants.extend(
                                    reactant_mol.GetSubstructMatches(c_n_pattern)
                                )
                                c_o_matches_reactants.extend(
                                    reactant_mol.GetSubstructMatches(c_o_pattern)
                                )

                    # If there are more C-N or C-O bonds in product than in reactants, bonds were formed
                    if len(c_n_matches_product) > len(c_n_matches_reactants) or len(
                        c_o_matches_product
                    ) > len(c_o_matches_reactants):
                        print("Heteroatom bond formation detected")
                        heteroatom_bond_formation_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return heteroatom_bond_formation_found
