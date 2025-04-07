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
    This function detects N-dibenzylation as a key protection/modification step.

    N-dibenzylation is a common protection strategy for amines where two benzyl
    groups are added to a primary amine or one benzyl group is added to a secondary
    amine, resulting in a tertiary amine with two benzyl groups.
    """
    has_dibenzylation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_dibenzylation

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is an N-alkylation reaction
                is_alkylation = checker.check_reaction(
                    "N-alkylation of primary amines with alkyl halides", rsmi
                ) or checker.check_reaction(
                    "N-alkylation of secondary amines with alkyl halides", rsmi
                )

                if is_alkylation:
                    # Check for primary or secondary amine in reactants
                    amine_reactant = None
                    for reactant in reactants_smiles:
                        if checker.check_fg(
                            "Primary amine", reactant
                        ) or checker.check_fg("Secondary amine", reactant):
                            amine_reactant = reactant
                            break

                    # Check for benzyl halides in reactants
                    benzyl_halides = []
                    for reactant in reactants_smiles:
                        if (
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                        ) and "c1ccccc1C" in reactant:
                            benzyl_halides.append(reactant)

                    # Check if product contains dibenzylamine pattern
                    product_has_dibenzyl = checker.check_fg(
                        "Tertiary amine", product_smiles
                    )

                    # Verify dibenzylation occurred
                    if (
                        amine_reactant
                        and len(benzyl_halides) >= 1
                        and product_has_dibenzyl
                    ):
                        # Count benzyl groups in product
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if product_mol:
                            benzyl_pattern = Chem.MolFromSmarts("c1ccccc1C[N]")
                            if product_mol.HasSubstructMatch(benzyl_pattern):
                                matches = product_mol.GetSubstructMatches(
                                    benzyl_pattern
                                )
                                if (
                                    len(matches) >= 2
                                ):  # At least two benzyl groups attached to nitrogen
                                    print(f"Detected N-dibenzylation at depth {depth}")
                                    print(f"Reaction SMILES: {rsmi}")
                                    has_dibenzylation = True

                # Also check for a single reaction that adds two benzyl groups at once
                if not has_dibenzylation:
                    for reactant in reactants_smiles:
                        if checker.check_fg("Primary amine", reactant):
                            # Check if product is a tertiary amine with two benzyl groups
                            product_mol = Chem.MolFromSmiles(product_smiles)
                            if product_mol and checker.check_fg(
                                "Tertiary amine", product_smiles
                            ):
                                benzyl_pattern = Chem.MolFromSmarts("c1ccccc1C[N]")
                                if product_mol.HasSubstructMatch(benzyl_pattern):
                                    matches = product_mol.GetSubstructMatches(
                                        benzyl_pattern
                                    )
                                    if len(matches) >= 2:
                                        print(
                                            f"Detected direct N-dibenzylation at depth {depth}"
                                        )
                                        print(f"Reaction SMILES: {rsmi}")
                                        has_dibenzylation = True
                                        break
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_dibenzylation
