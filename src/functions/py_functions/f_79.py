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
    This function detects late-stage amide formation from an amine and acid chloride.
    """
    # Initialize tracking variable
    has_amide_formation = False

    def dfs_traverse(node):
        nonlocal has_amide_formation

        if (
            node["type"] == "reaction" and not has_amide_formation
        ):  # Only check if not already found
            # Late stage corresponds to low depth
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                reactants = [
                    Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r
                ]
                product = Chem.MolFromSmiles(product_str) if product_str else None

                if product and any(reactants):
                    # Check for amide formation
                    amide_pattern = Chem.MolFromSmarts("[#6]-[#7]-[#6](=[#8])")
                    acid_chloride_pattern = Chem.MolFromSmarts("[#6](=[#8])[#17]")
                    amine_pattern = Chem.MolFromSmarts("[#6]-[#7]")

                    if (
                        product
                        and amide_pattern
                        and product.HasSubstructMatch(amide_pattern)
                    ):
                        has_acid_chloride = False
                        has_amine = False

                        for reactant in reactants:
                            if (
                                reactant
                                and acid_chloride_pattern
                                and reactant.HasSubstructMatch(acid_chloride_pattern)
                            ):
                                has_acid_chloride = True
                            if (
                                reactant
                                and amine_pattern
                                and reactant.HasSubstructMatch(amine_pattern)
                            ):
                                has_amine = True

                        if has_acid_chloride and has_amine:
                            has_amide_formation = True
                            print("Found late-stage amide formation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage amide formation detected: {has_amide_formation}")
    return has_amide_formation
