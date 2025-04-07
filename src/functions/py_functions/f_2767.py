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
    This function detects linear synthesis strategy with heterocycle formation
    followed by sequential functionalization.
    """
    # Track key transformations and their depths
    heterocycle_formation = {"detected": False, "depth": -1}
    functionalizations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a linear step (1-2 reactants)
                is_linear = len(reactants_smiles) <= 2

                if not is_linear:
                    return

                try:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product = Chem.MolFromSmiles(product_smiles)

                    if product is None or any(r is None for r in reactants):
                        return

                    # Check for heterocycle formation
                    heterocycle_patterns = [
                        Chem.MolFromSmarts("[n]1[n]cc[c]1"),  # pyrazole
                        Chem.MolFromSmarts("[o]1cc[n]c1"),  # oxazole
                        Chem.MolFromSmarts("[s]1cc[n]c1"),  # thiazole
                        Chem.MolFromSmarts("[n]1cc[n]c1"),  # imidazole
                        Chem.MolFromSmarts("[n]1ccc[n]1"),  # pyrazine
                        Chem.MolFromSmarts("[n]1cccc1"),  # pyrrole
                    ]

                    for pattern in heterocycle_patterns:
                        if product.HasSubstructMatch(pattern) and not any(
                            r.HasSubstructMatch(pattern)
                            for r in reactants
                            if r is not None
                        ):
                            heterocycle_formation["detected"] = True
                            heterocycle_formation["depth"] = depth
                            print(f"Detected heterocycle formation at depth {depth}")
                            break

                    # Check for functionalization reactions
                    functionalization_patterns = {
                        "halogenation": Chem.MolFromSmarts("[c,n]-[#9,#17,#35,#53]"),
                        "acylation": Chem.MolFromSmarts("[c,n,N]-[C](=[O])-[*]"),
                        "alkylation": Chem.MolFromSmarts("[c,n,N]-[C]"),
                        "oxidation": Chem.MolFromSmarts("[c,n,N]-[O]"),
                    }

                    for func_name, pattern in functionalization_patterns.items():
                        if product.HasSubstructMatch(pattern) and not all(
                            r.HasSubstructMatch(pattern)
                            for r in reactants
                            if r is not None
                        ):
                            functionalizations.append(
                                {"type": func_name, "depth": depth}
                            )
                            print(f"Detected {func_name} at depth {depth}")

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have heterocycle formation followed by at least one functionalization
    if heterocycle_formation["detected"] and functionalizations:
        # Check if functionalizations occur after heterocycle formation
        post_heterocycle_functionalizations = [
            f
            for f in functionalizations
            if f["depth"]
            < heterocycle_formation[
                "depth"
            ]  # Remember depth increases as we go earlier in synthesis
        ]

        if post_heterocycle_functionalizations:
            print("Detected linear heterocycle formation followed by functionalization")
            return True

    return False
