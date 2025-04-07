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
    This function detects if the synthesis route contains a sequence of
    transformations from alcohol to chloride to amine.
    """
    # Track reactions in sequence
    reactions_sequence = []

    def dfs_traverse(node, depth=0):
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Store reaction info with depth
            reactions_sequence.append(
                {
                    "depth": depth,
                    "reactants": reactants,
                    "product": product,
                    "rsmi": rsmi,
                }
            )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth (ascending order - from late to early in synthesis)
    reactions_sequence.sort(key=lambda x: x["depth"])

    # Check for the sequence: alcohol → chloride → amine
    has_alcohol_to_chloride = False
    has_chloride_to_amine = False

    for i in range(len(reactions_sequence) - 1):
        current_rxn = reactions_sequence[i]
        next_rxn = reactions_sequence[i + 1]

        # Check alcohol to chloride
        current_product_mol = Chem.MolFromSmiles(current_rxn["product"])
        next_reactants_mols = [
            Chem.MolFromSmiles(r) for r in next_rxn["reactants"] if r
        ]

        if current_product_mol:
            chloride_pattern = Chem.MolFromSmarts("[#6][Cl]")
            if current_product_mol.HasSubstructMatch(chloride_pattern):
                # Check if previous step had alcohol
                alcohol_pattern = Chem.MolFromSmarts("[#6][OH]")
                if any(
                    mol and mol.HasSubstructMatch(alcohol_pattern)
                    for mol in [
                        Chem.MolFromSmiles(r) for r in current_rxn["reactants"] if r
                    ]
                ):
                    has_alcohol_to_chloride = True
                    print(
                        f"Alcohol to chloride conversion detected at depth {current_rxn['depth']}"
                    )

            # Check chloride to amine
            if has_alcohol_to_chloride:
                amine_pattern = Chem.MolFromSmarts("[#6][N]")
                if next_rxn["product"] and Chem.MolFromSmiles(
                    next_rxn["product"]
                ).HasSubstructMatch(amine_pattern):
                    # Check if current product has chloride that's used in next reaction
                    if any(
                        mol and mol.HasSubstructMatch(chloride_pattern)
                        for mol in next_reactants_mols
                    ):
                        has_chloride_to_amine = True
                        print(
                            f"Chloride to amine conversion detected at depth {next_rxn['depth']}"
                        )

    sequence_found = has_alcohol_to_chloride and has_chloride_to_amine
    if sequence_found:
        print("Alcohol → Chloride → Amine sequence detected")

    return sequence_found
