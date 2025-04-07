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
    Detects the masked amine strategy using azide as an intermediate:
    alcohol → sulfonate → azide → amine → amide

    In retrosynthetic analysis, we'll see this in reverse:
    amide → amine → azide → sulfonate → alcohol
    """
    # Track if we've found each step in the sequence
    found_alcohol_activation = False
    found_azide_formation = False
    found_azide_reduction = False
    found_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_alcohol_activation, found_azide_formation, found_azide_reduction, found_amide_coupling

        if node["type"] == "reaction":
            # Extract reactants and product - in retrosynthesis, product is what we start with
            try:
                rsmi = node["metadata"]["rsmi"]
                # In retrosynthesis: product (right side) -> reactants (left side)
                reactants_part = rsmi.split(">")[0]  # What we're transforming to
                product_part = rsmi.split(">")[
                    -1
                ]  # What we're breaking down (starting with)

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for amide coupling (amide → amine in retrosynthesis)
                if (
                    checker.check_fg("Primary amide", product_part)
                    or checker.check_fg("Secondary amide", product_part)
                    or checker.check_fg("Tertiary amide", product_part)
                ) and checker.check_fg("Primary amine", reactants_part):
                    print(f"Found amide coupling step (retro) at depth {depth}")
                    found_amide_coupling = True

                # Check for amide formation reactions
                if (
                    checker.check_reaction(
                        "Carboxylic acid with primary amine to amide", rsmi
                    )
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        rsmi,
                    )
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                    or checker.check_reaction("Acylation of primary amines", rsmi)
                    or checker.check_reaction("Acylation of secondary amines", rsmi)
                ):
                    print(f"Found amide formation reaction at depth {depth}")
                    found_amide_coupling = True

                # Check for azide reduction (amine → azide in retrosynthesis)
                if checker.check_fg("Primary amine", product_part) and checker.check_fg(
                    "Azide", reactants_part
                ):
                    print(f"Found azide reduction step (retro) at depth {depth}")
                    found_azide_reduction = True

                # Check for azide reduction reactions
                if checker.check_reaction(
                    "Azide to amine reduction (Staudinger)", rsmi
                ) or checker.check_reaction(
                    "Reduction of nitro groups to amines", rsmi
                ):
                    print(f"Found azide reduction reaction at depth {depth}")
                    found_azide_reduction = True

                # Check for azide formation (azide → sulfonate in retrosynthesis)
                if checker.check_fg("Azide", product_part) and (
                    checker.check_fg("Triflate", reactants_part)
                    or checker.check_fg("Mesylate", reactants_part)
                    or checker.check_fg("Tosylate", reactants_part)
                    or checker.check_fg("Primary halide", reactants_part)
                    or checker.check_fg("Secondary halide", reactants_part)
                    or checker.check_fg("Tertiary halide", reactants_part)
                ):
                    print(f"Found azide formation step (retro) at depth {depth}")
                    found_azide_formation = True

                # Check for azide formation reactions
                if (
                    checker.check_reaction("Formation of Azides from halogens", rsmi)
                    or checker.check_reaction(
                        "Formation of Azides from boronic acids", rsmi
                    )
                    or checker.check_reaction("Alcohol to azide", rsmi)
                    or checker.check_reaction("Amine to azide", rsmi)
                ):
                    print(f"Found azide formation reaction at depth {depth}")
                    found_azide_formation = True

                # Check for alcohol activation (sulfonate → alcohol in retrosynthesis)
                if (
                    checker.check_fg("Triflate", product_part)
                    or checker.check_fg("Mesylate", product_part)
                    or checker.check_fg("Tosylate", product_part)
                ) and (
                    checker.check_fg("Primary alcohol", reactants_part)
                    or checker.check_fg("Secondary alcohol", reactants_part)
                    or checker.check_fg("Tertiary alcohol", reactants_part)
                    or checker.check_fg("Aromatic alcohol", reactants_part)
                ):
                    print(f"Found alcohol activation step (retro) at depth {depth}")
                    found_alcohol_activation = True

                # Check for alcohol to sulfonate reactions
                if (
                    checker.check_reaction("Formation of Sulfonic Esters", rsmi)
                    or checker.check_reaction("Alcohol to triflate conversion", rsmi)
                    or checker.check_reaction(
                        "Formation of Sulfonic Esters on TMS protected alcohol", rsmi
                    )
                    or checker.check_reaction(
                        "Alcohol to chloride_sulfonyl chloride", rsmi
                    )
                ):
                    print(f"Found alcohol to sulfonate reaction at depth {depth}")
                    found_alcohol_activation = True

                # Check for direct alcohol to halide conversions (which can be used instead of sulfonates)
                if (
                    checker.check_reaction("Alcohol to chloride_SOCl2", rsmi)
                    or checker.check_reaction("Alcohol to chloride_CHCl3", rsmi)
                    or checker.check_reaction("Alcohol to chloride_CH2Cl2", rsmi)
                    or checker.check_reaction("Alcohol to chloride_PCl5_ortho", rsmi)
                    or checker.check_reaction("Alcohol to chloride_POCl3_ortho", rsmi)
                    or checker.check_reaction("Alcohol to chloride_POCl3_para", rsmi)
                    or checker.check_reaction("Alcohol to chloride_POCl3", rsmi)
                    or checker.check_reaction("Alcohol to chloride_HCl", rsmi)
                    or checker.check_reaction("Alcohol to chloride_Salt", rsmi)
                    or checker.check_reaction("Alcohol to chloride_Other", rsmi)
                    or checker.check_reaction("Appel reaction", rsmi)
                ):
                    print(f"Found alcohol to halide conversion at depth {depth}")
                    found_alcohol_activation = True

                # Check for direct transformations that might skip steps
                # For example, direct alcohol to amine conversion
                if (
                    checker.check_fg("Primary alcohol", product_part)
                    or checker.check_fg("Secondary alcohol", product_part)
                    or checker.check_fg("Tertiary alcohol", product_part)
                ) and checker.check_fg("Primary amine", reactants_part):
                    print(f"Found direct alcohol to amine conversion at depth {depth}")
                    # This could potentially count as multiple steps
                    found_alcohol_activation = True
                    found_azide_formation = True
                    found_azide_reduction = True

                # Check for direct alcohol to azide conversion
                if (
                    checker.check_fg("Primary alcohol", product_part)
                    or checker.check_fg("Secondary alcohol", product_part)
                    or checker.check_fg("Tertiary alcohol", product_part)
                ) and checker.check_fg("Azide", reactants_part):
                    print(f"Found direct alcohol to azide conversion at depth {depth}")
                    found_alcohol_activation = True
                    found_azide_formation = True

                # Check for direct halide to azide conversion
                if (
                    checker.check_fg("Primary halide", product_part)
                    or checker.check_fg("Secondary halide", product_part)
                    or checker.check_fg("Tertiary halide", product_part)
                ) and checker.check_fg("Azide", reactants_part):
                    print(f"Found direct halide to azide conversion at depth {depth}")
                    found_azide_formation = True

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we found the complete sequence
    if (
        found_alcohol_activation
        and found_azide_formation
        and found_azide_reduction
        and found_amide_coupling
    ):
        print("Detected complete masked amine via azide strategy")
        return True
    else:
        print(
            f"Incomplete masked amine strategy: activation={found_alcohol_activation}, azide={found_azide_formation}, reduction={found_azide_reduction}, amide={found_amide_coupling}"
        )
        return False
