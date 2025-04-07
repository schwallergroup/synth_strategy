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
    Detects if the synthesis route uses a late-stage amide formation as the final step.
    """
    final_product = None
    amide_formation_detected = False

    # First identify the final product (root of the tree)
    if route["type"] == "mol":
        final_product = route["smiles"]
        print(f"Final product identified: {final_product}")

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_detected

        # Skip further processing if we already found an amide formation
        if amide_formation_detected:
            return

        # Check if this is a reaction node at a late stage (depth <= 2)
        if node["type"] == "reaction" and depth <= 2:
            print(f"Examining reaction step at depth {depth}")

            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Reaction SMILES: {rsmi}")
                print(f"Reactants: {reactants_smiles}")
                print(f"Product: {product_smiles}")

                # Check if this is an amide formation reaction using the checker function
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with ammonia to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with primary amine to imide",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with ammonia to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                ]

                is_amide_formation = any(
                    checker.check_reaction(rxn_type, rsmi)
                    for rxn_type in amide_formation_reactions
                )

                if is_amide_formation:
                    print(f"Detected amide formation reaction at depth {depth}")
                    # We've confirmed this is an amide formation reaction
                    amide_formation_detected = True
                    return

                # If reaction type check failed, fall back to functional group analysis
                if not is_amide_formation:
                    # Parse molecules to identify the specific atoms involved
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Check for amine in reactants
                    amine_reactant_idx = None
                    for idx, r in enumerate(reactants_smiles):
                        if checker.check_fg("Primary amine", r) or checker.check_fg(
                            "Secondary amine", r
                        ):
                            amine_reactant_idx = idx
                            break

                    # Check for carbonyl-containing reactant
                    carbonyl_reactant_idx = None
                    for idx, r in enumerate(reactants_smiles):
                        if (
                            checker.check_fg("Carboxylic acid", r)
                            or checker.check_fg("Acyl halide", r)
                            or checker.check_fg("Ester", r)
                            or checker.check_fg("Anhydride", r)
                        ):
                            carbonyl_reactant_idx = idx
                            break

                    # Check if both reactants are present
                    if (
                        amine_reactant_idx is not None
                        and carbonyl_reactant_idx is not None
                    ):
                        # Check if product has a new amide bond
                        has_amide = (
                            checker.check_fg("Primary amide", product_smiles)
                            or checker.check_fg("Secondary amide", product_smiles)
                            or checker.check_fg("Tertiary amide", product_smiles)
                        )

                        # Check if the reaction involves acylation of an amine
                        # Look for patterns where an amine reacts with a carbonyl compound
                        if has_amide:
                            # Check if the amine is primary or secondary
                            amine_reactant = reactants_smiles[amine_reactant_idx]
                            has_primary_amine = checker.check_fg(
                                "Primary amine", amine_reactant
                            )
                            has_secondary_amine = checker.check_fg(
                                "Secondary amine", amine_reactant
                            )

                            # Check if the carbonyl compound is an acyl halide, acid, ester, etc.
                            carbonyl_reactant = reactants_smiles[carbonyl_reactant_idx]
                            has_acyl_halide = checker.check_fg(
                                "Acyl halide", carbonyl_reactant
                            )
                            has_carboxylic_acid = checker.check_fg(
                                "Carboxylic acid", carbonyl_reactant
                            )
                            has_ester = checker.check_fg("Ester", carbonyl_reactant)

                            print(
                                f"Amine reactant: {has_primary_amine or has_secondary_amine}"
                            )
                            print(
                                f"Carbonyl reactant: {has_acyl_halide or has_carboxylic_acid or has_ester}"
                            )

                            # If we have the right reactants and product has an amide, it's likely an amide formation
                            if (has_primary_amine or has_secondary_amine) and (
                                has_acyl_halide or has_carboxylic_acid or has_ester
                            ):
                                print(
                                    f"Detected amide formation at depth {depth} through functional group analysis"
                                )
                                amide_formation_detected = True
                                return
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    print(f"Final result: {amide_formation_detected}")
    return amide_formation_detected
