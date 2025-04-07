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
    Detects if the synthesis route involves sequential elaboration of fragments
    before final coupling.
    """
    # Track fragment modifications and couplings
    fragment_modifications = {}  # Track modifications by fragment SMILES
    fragment_lineage = {}  # Track product -> reactant relationships
    coupling_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is a coupling reaction (2+ significant reactants)
            if len(reactants_smiles) >= 2:
                significant_fragments = 0
                for smi in reactants_smiles:
                    mol = Chem.MolFromSmiles(smi)
                    if mol and (
                        Chem.GetSSSR(mol) or mol.GetNumHeavyAtoms() > 5
                    ):  # Has at least one ring or is substantial
                        significant_fragments += 1

                if significant_fragments >= 2:
                    coupling_reactions.append((depth, rsmi, reactants_smiles))

                    # Track which fragments were used in this coupling
                    for reactant in reactants_smiles:
                        if reactant in fragment_modifications:
                            fragment_modifications[reactant].append(
                                (depth, rsmi, "coupled")
                            )

            # For all reactions, track reactant -> product relationships
            for reactant in reactants_smiles:
                if reactant not in fragment_modifications:
                    fragment_modifications[reactant] = []
                fragment_modifications[reactant].append((depth, rsmi, "modified"))

                # Track lineage
                if product_smiles not in fragment_lineage:
                    fragment_lineage[product_smiles] = []
                fragment_lineage[product_smiles].append(reactant)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Analyze if modifications occur before couplings
    has_sequential_elaboration = False

    if coupling_reactions:
        # Get the depth of the first coupling reaction
        first_coupling_depth = min(depth for depth, _, _ in coupling_reactions)

        # Check if there are modifications at greater depths (earlier in synthesis)
        deeper_modifications = [
            depth
            for frag, mods in fragment_modifications.items()
            for depth, _, mod_type in mods
            if depth > first_coupling_depth and mod_type == "modified"
        ]

        # We have sequential elaboration if modifications happen before coupling
        has_sequential_elaboration = len(deeper_modifications) > 0

        print(f"First coupling depth: {first_coupling_depth}")
        print(f"Deeper modifications: {len(deeper_modifications)}")

    print(f"Sequential fragment elaboration strategy: {has_sequential_elaboration}")
    print(
        f"Coupling reactions: {len(coupling_reactions)}, Fragment modifications: {sum(len(mods) for mods in fragment_modifications.values())}"
    )

    return has_sequential_elaboration
