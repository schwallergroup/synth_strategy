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
    This function detects a synthetic strategy involving N-methylation reactions
    in the late stages of synthesis (low depth in retrosynthetic tree).
    """
    # Track N-methylations at different depths
    n_methylations_by_depth = {}
    max_depth = 0

    def is_n_methylation(rsmi):
        """Check if a reaction is an N-methylation using checker functions"""
        # Check for various N-methylation reaction types
        n_methylation_reactions = [
            "N-methylation",
            "Eschweiler-Clarke Primary Amine Methylation",
            "Eschweiler-Clarke Secondary Amine Methylation",
            "Reductive methylation of primary amine with formaldehyde",
            "Methylation with MeI_primary",
            "Methylation with MeI_secondary",
            "Methylation with MeI_tertiary",
            "DMS Amine methylation",
            "Parnes methylation",
            "Methylation with DMS",
            "Methylation",  # Check general methylation first, then verify if it's on nitrogen
        ]

        for reaction_type in n_methylation_reactions:
            if checker.check_reaction(reaction_type, rsmi):
                # For general methylation, verify it's specifically N-methylation
                if reaction_type == "Methylation":
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Find the main reactant (usually the largest molecule)
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                        product_mol = Chem.MolFromSmiles(product)

                        if not product_mol or not all(reactant_mols):
                            continue

                        main_reactant = max(
                            reactant_mols, key=lambda m: m.GetNumAtoms()
                        )

                        # Check if the product has more N-CH3 groups than the main reactant
                        # Use a more general pattern for N-CH3 groups
                        n_methyl_pattern = Chem.MolFromSmarts("[#7][CH3]")

                        reactant_n_methyl_count = len(
                            main_reactant.GetSubstructMatches(n_methyl_pattern)
                        )
                        product_n_methyl_count = len(
                            product_mol.GetSubstructMatches(n_methyl_pattern)
                        )

                        # Check if there's an increase in N-CH3 groups
                        if product_n_methyl_count > reactant_n_methyl_count:
                            # Verify that the reactant has nitrogen atoms that could be methylated
                            nitrogen_pattern = Chem.MolFromSmarts("[#7]")
                            if main_reactant.HasSubstructMatch(nitrogen_pattern):
                                print(
                                    f"Detected general methylation on nitrogen: {rsmi}"
                                )
                                return True
                    except Exception as e:
                        print(f"Error in methylation verification: {e}")
                        continue
                else:
                    # For specific N-methylation reactions, we trust the checker
                    print(f"Detected {reaction_type} reaction: {rsmi}")
                    return True

        # Check for functional groups that might indicate N-methylation
        try:
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant has primary/secondary amine and product has tertiary amine or N-CH3
            for reactant in reactants:
                if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                    "Secondary amine", reactant
                ):
                    if checker.check_fg(
                        "Tertiary amine", product
                    ) or Chem.MolFromSmiles(product).HasSubstructMatch(
                        Chem.MolFromSmarts("[#7][CH3]")
                    ):
                        # Verify increase in N-CH3 groups
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        product_mol = Chem.MolFromSmiles(product)
                        n_methyl_pattern = Chem.MolFromSmarts("[#7][CH3]")

                        if reactant_mol and product_mol:
                            reactant_n_methyl_count = len(
                                reactant_mol.GetSubstructMatches(n_methyl_pattern)
                            )
                            product_n_methyl_count = len(
                                product_mol.GetSubstructMatches(n_methyl_pattern)
                            )

                            if product_n_methyl_count > reactant_n_methyl_count:
                                print(
                                    f"Detected N-methylation via functional group analysis: {rsmi}"
                                )
                                return True
        except Exception as e:
            print(f"Error in functional group analysis: {e}")

        return False

    def dfs_traverse(node, depth=0):
        nonlocal n_methylations_by_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            try:
                # Extract reaction SMILES
                rsmi = node["metadata"]["rsmi"]

                # Check if this is an N-methylation reaction
                if is_n_methylation(rsmi):
                    if depth not in n_methylations_by_depth:
                        n_methylations_by_depth[depth] = 0
                    n_methylations_by_depth[depth] += 1
                    print(f"Found N-methylation at depth {depth}: {rsmi}")
            except KeyError:
                print(f"Missing rsmi in reaction metadata")
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting route traversal to find N-methylation reactions")
    dfs_traverse(route)

    # Calculate the percentage of N-methylations in the first half of the synthesis
    total_methylations = sum(n_methylations_by_depth.values())
    if total_methylations == 0:
        print("No N-methylation reactions found in the route")
        return False

    # Define "late stage" as the first half of the synthesis (lower depths)
    late_stage_threshold = max_depth / 2
    late_stage_methylations = sum(
        count
        for depth, count in n_methylations_by_depth.items()
        if depth <= late_stage_threshold
    )

    # Strategy criteria: at least 2 N-methylations with at least 50% in late stage
    late_stage_percentage = (
        late_stage_methylations / total_methylations if total_methylations > 0 else 0
    )
    strategy_detected = total_methylations >= 2 and late_stage_percentage >= 0.5

    print(f"Late-stage N-methylation strategy detected: {strategy_detected}")
    print(f"Total N-methylations: {total_methylations}")
    print(
        f"Late-stage N-methylations: {late_stage_methylations} ({late_stage_percentage:.1%})"
    )
    print(f"Max depth: {max_depth}, Late-stage threshold: {late_stage_threshold}")
    print(f"N-methylations by depth: {n_methylations_by_depth}")

    return strategy_detected
