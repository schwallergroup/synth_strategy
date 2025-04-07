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
    Detects a strategy involving late-stage C-N coupling between an aryl halide and
    a nitrogen heterocycle.
    """
    c_n_coupling_detected = False

    # List of nitrogen heterocycles to check
    n_heterocycles = [
        "pyridine",
        "pyrrole",
        "imidazole",
        "pyrazole",
        "triazole",
        "tetrazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "indole",
        "quinoline",
        "isoquinoline",
        "benzimidazole",
        "benzotriazole",
    ]

    # List of C-N coupling reactions to check
    cn_coupling_reactions = [
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "N-arylation_heterocycles",
        "Goldberg coupling",
        "Ullmann-Goldberg Substitution amine",
        "Goldberg coupling aryl amine-aryl chloride",
        "Goldberg coupling aryl amide-aryl chloride",
        "Ullmann condensation",
    ]

    # Define late stage depth threshold
    LATE_STAGE_DEPTH = 2  # Consider depths 0, 1, and 2 as late stage

    def dfs_traverse(node, depth=0):
        nonlocal c_n_coupling_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Extract depth from metadata if available
            node_depth = depth
            if "ID" in node.get("metadata", {}):
                depth_match = re.search(r"Depth: (\d+)", node.get("metadata", {}).get("ID", ""))
                if depth_match:
                    node_depth = int(depth_match.group(1))

            # Check if this is a late-stage reaction
            if node_depth <= LATE_STAGE_DEPTH:
                print(f"Examining late-stage reaction at depth {node_depth}: {rsmi}")

                # Find reactants with aryl halides
                aryl_halide_reactants = [
                    r for r in reactants if r and checker.check_fg("Aromatic halide", r)
                ]

                # Find reactants with nitrogen heterocycles
                n_heterocycle_reactants = [
                    r
                    for r in reactants
                    if r and any(checker.check_ring(ring, r) for ring in n_heterocycles)
                ]

                # Check if we have both required components
                has_aryl_halide = len(aryl_halide_reactants) > 0
                has_n_heterocycle = len(n_heterocycle_reactants) > 0

                print(
                    f"Aryl halide reactants: {has_aryl_halide}, N-heterocycle reactants: {has_n_heterocycle}"
                )

                # First check if this is a known C-N coupling reaction
                is_cn_coupling = any(
                    checker.check_reaction(rxn, rsmi) for rxn in cn_coupling_reactions
                )
                print(f"Is known C-N coupling reaction: {is_cn_coupling}")

                # If not a known reaction type, check for the chemical transformation directly
                if not is_cn_coupling and has_aryl_halide and has_n_heterocycle:
                    # Check if the product contains copper - might be a copper-mediated coupling
                    copper_present = "[Cu]" in rsmi

                    # Check if the product contains both the aryl group and the N-heterocycle
                    product_has_n_heterocycle = any(
                        checker.check_ring(ring, product) for ring in n_heterocycles
                    )

                    # Check if the aryl halide is consumed (no longer present in product)
                    product_has_aryl_halide = checker.check_fg("Aromatic halide", product)

                    # If copper is present, N-heterocycle is in product, and aryl halide is consumed,
                    # this is likely a C-N coupling
                    if copper_present and product_has_n_heterocycle and not product_has_aryl_halide:
                        print(
                            f"Detected copper-mediated C-N coupling based on reactants and products"
                        )
                        is_cn_coupling = True

                # If we have all components and it's a C-N coupling reaction
                if has_aryl_halide and has_n_heterocycle and is_cn_coupling:
                    # Verify that the product contains both the aryl group and the N-heterocycle
                    product_has_n_heterocycle = any(
                        checker.check_ring(ring, product) for ring in n_heterocycles
                    )

                    # For a successful C-N coupling, the product should no longer have an aryl halide
                    # but should still have an aromatic ring
                    product_has_aromatic = (
                        checker.check_fg("Phenol", product)
                        or checker.check_fg("Aniline", product)
                        or checker.check_ring("benzene", product)
                        or any(
                            checker.check_ring(ring, product)
                            for ring in ["naphthalene", "anthracene"]
                        )
                    )

                    print(f"Product has N-heterocycle: {product_has_n_heterocycle}")
                    print(f"Product has aromatic ring: {product_has_aromatic}")

                    if product_has_n_heterocycle and product_has_aromatic:
                        print(f"Detected C-N coupling reaction at late stage (depth {node_depth})")
                        c_n_coupling_detected = True

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return c_n_coupling_detected
