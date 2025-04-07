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
    This function detects if the synthetic route follows a linear fragment assembly strategy
    with sequential heteroatom substitutions.
    """
    # Track the number of sequential heteroatom substitutions
    heteroatom_substitutions = 0
    # Track if the route is linear (each reaction has only one product)
    is_linear = True

    def dfs_traverse(node):
        nonlocal heteroatom_substitutions, is_linear

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Check if there's only one product (linear synthesis)
                if "." in product_part:
                    is_linear = False

                try:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                    product = Chem.MolFromSmiles(product_part)

                    if all(r is not None for r in reactants) and product is not None:
                        # Check for C-O or C-N bond formation
                        c_o_pattern = Chem.MolFromSmarts("[#6]-[#8]")
                        c_n_pattern = Chem.MolFromSmarts("[#6]-[#7]")

                        # Count C-O and C-N bonds in reactants and product
                        reactants_c_o_count = sum(
                            len(r.GetSubstructMatches(c_o_pattern)) for r in reactants
                        )
                        reactants_c_n_count = sum(
                            len(r.GetSubstructMatches(c_n_pattern)) for r in reactants
                        )

                        product_c_o_count = len(product.GetSubstructMatches(c_o_pattern))
                        product_c_n_count = len(product.GetSubstructMatches(c_n_pattern))

                        # Check if C-O or C-N bonds increased
                        if (
                            product_c_o_count > reactants_c_o_count
                            or product_c_n_count > reactants_c_n_count
                        ):
                            heteroatom_substitutions += 1
                except:
                    print(f"Error processing reaction SMILES: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if the route is linear and has multiple heteroatom substitutions
    return is_linear and heteroatom_substitutions >= 2
