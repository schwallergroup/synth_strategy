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
    This function detects a linear synthesis strategy where each step builds on a single precursor.
    It checks if most reactions have only one product-contributing reactant by analyzing structural similarity.
    """
    reaction_count = 0
    linear_reaction_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, linear_reaction_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            reaction_count += 1
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Handle multiple reactants
            reactants = reactants_part.split(".")
            product_mol = Chem.MolFromSmiles(product_part)

            if not product_mol:
                print(f"Warning: Could not parse product SMILES at depth {depth}: {product_part}")
                return

            # Find the reactant with the highest structural similarity to the product
            max_similarity = 0
            main_contributor = None

            for reactant in reactants:
                try:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    # Skip very small molecules (likely reagents, not main scaffold)
                    if reactant_mol.GetNumAtoms() < 3:
                        continue

                    # Calculate structural similarity using Maximum Common Substructure
                    mcs = rdFMCS.FindMCS(
                        [reactant_mol, product_mol],
                        atomCompare=rdFMCS.AtomCompare.CompareElements,
                        bondCompare=rdFMCS.BondCompare.CompareOrder,
                        completeRingsOnly=True,
                        timeout=1,
                    )

                    if mcs.numAtoms > 0:
                        # Calculate similarity as proportion of product atoms in the MCS
                        similarity = mcs.numAtoms / product_mol.GetNumAtoms()

                        if similarity > max_similarity:
                            max_similarity = similarity
                            main_contributor = reactant
                except Exception as e:
                    print(f"Error processing reactant at depth {depth}: {e}")

            # If one reactant contributes significantly to the product structure (>50% similarity)
            # consider it a linear step
            if max_similarity > 0.5:
                linear_reaction_count += 1
                print(
                    f"Linear reaction step found at depth {depth} (similarity: {max_similarity:.2f})"
                )
            else:
                print(
                    f"Non-linear reaction step at depth {depth} (max similarity: {max_similarity:.2f})"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if most reactions (>75%) are linear
    linear_ratio = linear_reaction_count / reaction_count if reaction_count > 0 else 0
    print(f"Linear reactions: {linear_reaction_count}/{reaction_count} ({linear_ratio:.2f})")

    return reaction_count > 0 and linear_ratio > 0.75
