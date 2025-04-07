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
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    This function detects if the synthesis route involves early-stage heterocycle formation,
    specifically looking for nitrogen-containing heterocycle construction.
    """
    heterocycle_formation_detected = False
    max_depth_with_heterocycle = -1

    # List of nitrogen-containing heterocycles to check
    heterocycles = [
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_detected, max_depth_with_heterocycle

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = (
                    Chem.MolFromSmiles(product_smiles) if product_smiles else None
                )

                if all(reactant_mols) and product_mol:
                    # Count rings in reactants and product
                    reactant_ring_counts = [
                        mol.GetRingInfo().NumRings() for mol in reactant_mols
                    ]
                    product_ring_count = product_mol.GetRingInfo().NumRings()

                    # Check if product has more rings than any reactant
                    if product_ring_count > max(reactant_ring_counts, default=0):
                        # Check if any of the specified heterocycles are in the product
                        product_heterocycles = [
                            ring
                            for ring in heterocycles
                            if checker.check_ring(ring, product_smiles)
                        ]

                        if product_heterocycles:
                            # Check if these heterocycles are not present in any reactant
                            new_heterocycle_formed = True
                            for ring in product_heterocycles:
                                if any(
                                    checker.check_ring(ring, r)
                                    for r in reactants_smiles
                                ):
                                    new_heterocycle_formed = False
                                    break

                            if new_heterocycle_formed:
                                heterocycle_formation_detected = True
                                max_depth_with_heterocycle = max(
                                    max_depth_with_heterocycle, depth
                                )
                                print(
                                    f"Detected heterocycle formation at depth {depth}: {', '.join(product_heterocycles)}"
                                )
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if heterocycle formation occurred in early stage (high depth)
    early_stage = (
        max_depth_with_heterocycle >= 4
    )  # Considering depth ≥ 4 as early stage

    return heterocycle_formation_detected and early_stage
