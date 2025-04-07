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
    This function detects a strategy involving late-stage amide coupling
    as the final step in the synthesis.
    """
    has_late_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_amide_coupling

        # Check if this is a reaction node at a late stage (depth <= 2)
        if node["type"] == "reaction" and depth <= 2:
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                reagents_smiles = (
                    rsmi.split(">")[1].split(".") if len(rsmi.split(">")) > 2 else []
                )
                product_smiles = rsmi.split(">")[-1]

                print(
                    f"Checking potential late-stage reaction at depth {depth}: {rsmi}"
                )

                # Check if this is an amide coupling reaction using the checker functions
                is_amide_coupling = any(
                    [
                        checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                            rsmi,
                        ),
                        checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                            rsmi,
                        ),
                        checker.check_reaction(
                            "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                            rsmi,
                        ),
                        checker.check_reaction(
                            "Carboxylic acid with primary amine to amide", rsmi
                        ),
                        checker.check_reaction(
                            "Ester with primary amine to amide", rsmi
                        ),
                        checker.check_reaction(
                            "Ester with secondary amine to amide", rsmi
                        ),
                        checker.check_reaction(
                            "Acyl chloride with secondary amine to amide", rsmi
                        ),
                        checker.check_reaction(
                            "Acyl chloride with ammonia to amide", rsmi
                        ),
                        checker.check_reaction("Ester with ammonia to amide", rsmi),
                        checker.check_reaction("Schotten-Baumann_amide", rsmi),
                        checker.check_reaction("Acylation of primary amines", rsmi),
                        checker.check_reaction("Acylation of secondary amines", rsmi),
                        checker.check_reaction(
                            "Carboxylic acid to amide conversion", rsmi
                        ),
                    ]
                )

                # If not directly identified, check for common coupling patterns
                if not is_amide_coupling:
                    # Check for presence of coupling reagents in the reaction
                    coupling_reagents_present = any(
                        ["CCN=C=NCCCN(C)C" in reagent for reagent in reagents_smiles]
                    ) or any(
                        ["CCN=C=NCCCN(C)C" in reactant for reactant in reactants_smiles]
                    )

                    # Check for amide formation pattern manually
                    has_carboxylic_acid = any(
                        [
                            checker.check_fg("Carboxylic acid", r)
                            for r in reactants_smiles
                        ]
                    )
                    has_amine = any(
                        [
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            or checker.check_fg("Tertiary amine", r)
                            for r in reactants_smiles
                        ]
                    )

                    # Check if product has an amide
                    has_amide_in_product = any(
                        [
                            checker.check_fg("Secondary amide", product_smiles),
                            checker.check_fg("Primary amide", product_smiles),
                            checker.check_fg("Tertiary amide", product_smiles),
                        ]
                    )

                    # If we have coupling reagents, acid, amine, and amide in product, it's likely an amide coupling
                    if (
                        coupling_reagents_present
                        and has_carboxylic_acid
                        and has_amine
                        and has_amide_in_product
                    ):
                        is_amide_coupling = True
                        print(
                            f"Detected amide coupling with coupling reagents at depth {depth}"
                        )

                if is_amide_coupling:
                    print(f"Found amide coupling reaction at depth {depth}")

                    # Verify that an amide is formed in the product
                    has_amide_in_product = any(
                        [
                            checker.check_fg("Secondary amide", product_smiles),
                            checker.check_fg("Primary amide", product_smiles),
                            checker.check_fg("Tertiary amide", product_smiles),
                        ]
                    )

                    if has_amide_in_product:
                        print(f"Confirmed amide in product: {product_smiles}")

                        # Check if reactants have the necessary components for amide coupling
                        has_amine = False
                        has_carboxylic_acid_or_derivative = False

                        for reactant in reactants_smiles:
                            if (
                                checker.check_fg("Primary amine", reactant)
                                or checker.check_fg("Secondary amine", reactant)
                                or checker.check_fg("Tertiary amine", reactant)
                                or checker.check_fg("Aniline", reactant)
                            ):
                                has_amine = True
                                print(f"Found amine in reactant: {reactant}")

                            if (
                                checker.check_fg("Carboxylic acid", reactant)
                                or checker.check_fg("Acyl halide", reactant)
                                or checker.check_fg("Ester", reactant)
                                or checker.check_fg("Anhydride", reactant)
                            ):
                                has_carboxylic_acid_or_derivative = True
                                print(
                                    f"Found carboxylic acid or derivative in reactant: {reactant}"
                                )

                        # If we don't find both components but have coupling reagents, it might still be an amide coupling
                        if (
                            has_amine and has_carboxylic_acid_or_derivative
                        ) or coupling_reagents_present:
                            has_late_amide_coupling = True
                            print(
                                f"Confirmed late-stage amide coupling at depth {depth}"
                            )
                        else:
                            print(f"Missing required reactants for amide coupling")
                    else:
                        print(f"No amide found in product despite reaction match")
                else:
                    print(f"Not an amide coupling reaction at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: has_late_amide_coupling = {has_late_amide_coupling}")

    return has_late_amide_coupling
