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
    This function detects a synthesis strategy that involves halide displacement
    reactions, particularly the conversion of alkyl halides to other functional groups.
    """
    # Track halide displacement reactions
    halide_displacements = 0

    def dfs_traverse(node, depth=0):
        nonlocal halide_displacements

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
            product = Chem.MolFromSmiles(product_part)

            if product is None or any(r is None for r in reactants):
                print("Warning: Could not parse some molecules in reaction")
                return

            # Check for halide displacement patterns
            # Alkyl halide to nitrile
            if any(
                r.HasSubstructMatch(Chem.MolFromSmarts("[#6][Br,Cl,I,F]")) for r in reactants
            ) and product.HasSubstructMatch(Chem.MolFromSmarts("[#6]C#N")):
                halide_displacements += 1
                print(f"Detected halide displacement: Alkyl halide to nitrile at depth {depth}")

            # Alkyl halide to other nucleophilic substitution products
            if any(
                r.HasSubstructMatch(Chem.MolFromSmarts("[#6][Br,Cl,I,F]")) for r in reactants
            ) and (
                product.HasSubstructMatch(Chem.MolFromSmarts("[#6][O,N,S]"))
                or product.HasSubstructMatch(Chem.MolFromSmarts("[#6][#6]"))
            ):
                halide_displacements += 1
                print(
                    f"Detected halide displacement: General nucleophilic substitution at depth {depth}"
                )

        # Process children
        if "children" in node:
            for child in node["children"]:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Halide displacement reactions: {halide_displacements}")

    # Return True if we have at least one halide displacement reaction
    return halide_displacements >= 1
