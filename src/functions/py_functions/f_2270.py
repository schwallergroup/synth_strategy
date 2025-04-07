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
    This function detects a linear synthesis strategy where multiple functional group
    transformations occur on a side chain while maintaining a heterocyclic core structure.
    """
    # Track if we've found evidence of our strategy
    core_preserved = False
    side_chain_modifications = 0
    oxazoline_present = False

    # Define SMARTS patterns
    oxazoline_pattern = Chem.MolFromSmarts("[#6]1[#6][#7][#6][#6]1=N[O]")

    # Functional groups to track
    azide_pattern = Chem.MolFromSmarts("[#7-]=[#7+]=[#7]")
    tosylate_pattern = Chem.MolFromSmarts("[#8][S](=O)(=O)c1ccc(C)cc1")
    alcohol_pattern = Chem.MolFromSmarts("[#8H]")
    ester_pattern = Chem.MolFromSmarts("[#6][#8]C(=O)[#6]")

    # Track which functional groups we've seen
    seen_functional_groups = set()

    def dfs_traverse(node):
        nonlocal core_preserved, side_chain_modifications, oxazoline_present, seen_functional_groups

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for oxazoline core
                if mol.HasSubstructMatch(oxazoline_pattern):
                    oxazoline_present = True

                # Check for functional groups
                if mol.HasSubstructMatch(azide_pattern):
                    seen_functional_groups.add("azide")
                if mol.HasSubstructMatch(tosylate_pattern):
                    seen_functional_groups.add("tosylate")
                if mol.HasSubstructMatch(alcohol_pattern):
                    seen_functional_groups.add("alcohol")
                if mol.HasSubstructMatch(ester_pattern):
                    seen_functional_groups.add("ester")

        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if the core structure is preserved in this reaction
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and any(r for r in reactant_mols if r):
                reactant_has_core = any(
                    r.HasSubstructMatch(oxazoline_pattern) for r in reactant_mols if r
                )
                product_has_core = product_mol.HasSubstructMatch(oxazoline_pattern)

                if reactant_has_core and product_has_core:
                    core_preserved = True
                    side_chain_modifications += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Determine if our strategy is present
    strategy_present = (
        oxazoline_present
        and core_preserved
        and side_chain_modifications >= 3
        and len(seen_functional_groups) >= 3
    )

    print(f"Linear side chain modification strategy detected: {strategy_present}")
    print(f"Oxazoline present: {oxazoline_present}")
    print(f"Core preserved: {core_preserved}")
    print(f"Side chain modifications: {side_chain_modifications}")
    print(f"Functional groups seen: {seen_functional_groups}")

    return strategy_present
