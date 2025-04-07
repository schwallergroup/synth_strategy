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
    Detects a sequence of functional group interconversions leading to a key bond formation.
    Specifically looks for alcohol→halide→amine→amide transformation sequence.
    """
    # Track transformations with their depths
    transformations = {
        "alcohol_to_halide": None,
        "halide_to_amine": None,
        "amine_to_amide": None,
    }

    print(
        "Starting to analyze the synthesis route for functional group interconversion sequence"
    )

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for alcohol to halide transformation
                alcohol_in_reactants = any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    for r in reactants_smiles
                )

                halide_in_product = (
                    checker.check_fg("Primary halide", product_smiles)
                    or checker.check_fg("Secondary halide", product_smiles)
                    or checker.check_fg("Tertiary halide", product_smiles)
                )

                if alcohol_in_reactants and halide_in_product:
                    if (
                        checker.check_reaction("Alcohol to chloride_SOCl2", rsmi)
                        or checker.check_reaction(
                            "Alcohol to chloride_PCl5_ortho", rsmi
                        )
                        or checker.check_reaction("Alcohol to chloride_POCl3", rsmi)
                        or checker.check_reaction("Alcohol to chloride_HCl", rsmi)
                        or checker.check_reaction("Alcohol to chloride_Salt", rsmi)
                        or checker.check_reaction("Alcohol to chloride_Other", rsmi)
                        or checker.check_reaction("Alkyl bromides from alcohols", rsmi)
                        or checker.check_reaction("Alkyl iodides from alcohols", rsmi)
                        or checker.check_reaction("Appel reaction", rsmi)
                        or checker.check_reaction(
                            "Alcohol to chloride_sulfonyl chloride", rsmi
                        )
                        or checker.check_reaction(
                            "PBr3 and alcohol to alkyl bromide", rsmi
                        )
                    ):
                        transformations["alcohol_to_halide"] = depth
                        print(
                            f"Found alcohol to halide transformation at depth {depth}"
                        )
                    else:
                        # Fallback detection if no specific reaction type matches
                        transformations["alcohol_to_halide"] = depth
                        print(
                            f"Found alcohol to halide transformation at depth {depth} (fallback detection)"
                        )

                # Check for halide to amine transformation
                halide_in_reactants = any(
                    checker.check_fg("Primary halide", r)
                    or checker.check_fg("Secondary halide", r)
                    or checker.check_fg("Tertiary halide", r)
                    for r in reactants_smiles
                )

                amine_in_product = (
                    checker.check_fg("Primary amine", product_smiles)
                    or checker.check_fg("Secondary amine", product_smiles)
                    or checker.check_fg("Tertiary amine", product_smiles)
                )

                if halide_in_reactants and amine_in_product:
                    if (
                        checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction(
                            "N-alkylation of secondary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction("Alkylation of amines", rsmi)
                        or checker.check_reaction("N-methylation", rsmi)
                    ):
                        transformations["halide_to_amine"] = depth
                        print(f"Found halide to amine transformation at depth {depth}")
                    else:
                        # Fallback detection if no specific reaction type matches
                        transformations["halide_to_amine"] = depth
                        print(
                            f"Found halide to amine transformation at depth {depth} (fallback detection)"
                        )

                # Check for amine to amide transformation
                amine_in_reactants = any(
                    checker.check_fg("Primary amine", r)
                    or checker.check_fg("Secondary amine", r)
                    or checker.check_fg("Tertiary amine", r)
                    for r in reactants_smiles
                )

                amide_in_product = (
                    checker.check_fg("Primary amide", product_smiles)
                    or checker.check_fg("Secondary amide", product_smiles)
                    or checker.check_fg("Tertiary amide", product_smiles)
                )

                if amine_in_reactants and amide_in_product:
                    if (
                        checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                            rsmi,
                        )
                        or checker.check_reaction("Acylation of primary amines", rsmi)
                        or checker.check_reaction("Acylation of secondary amines", rsmi)
                        or checker.check_reaction(
                            "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                            rsmi,
                        )
                        or checker.check_reaction(
                            "Acyl chloride with secondary amine to amide", rsmi
                        )
                        or checker.check_reaction(
                            "Carboxylic acid with primary amine to amide", rsmi
                        )
                        or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                        or checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                            rsmi,
                        )
                    ):
                        transformations["amine_to_amide"] = depth
                        print(f"Found amine to amide transformation at depth {depth}")
                    else:
                        # Fallback detection if no specific reaction type matches
                        transformations["amine_to_amide"] = depth
                        print(
                            f"Found amine to amide transformation at depth {depth} (fallback detection)"
                        )

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Transformation depths: {transformations}")

    # Check if we observed the complete sequence in the correct order
    if all(v is not None for v in transformations.values()):
        # In retrosynthetic direction, higher depth means earlier in synthesis
        # So alcohol→halide→amine→amide means:
        # depth(alcohol_to_halide) > depth(halide_to_amine) > depth(amine_to_amide)
        if (
            transformations["alcohol_to_halide"]
            > transformations["halide_to_amine"]
            > transformations["amine_to_amide"]
        ):
            print(
                "Detected alcohol→halide→amine→amide transformation sequence in correct order"
            )
            return True
        else:
            print("All transformations found but not in the correct sequence order")
            print(
                f"Expected: alcohol_to_halide ({transformations['alcohol_to_halide']}) > "
                f"halide_to_amine ({transformations['halide_to_amine']}) > "
                f"amine_to_amide ({transformations['amine_to_amide']})"
            )
    else:
        print("Not all required transformations were found in the synthesis route")

    return False
