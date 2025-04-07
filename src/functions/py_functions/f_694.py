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
    This function detects late-stage amide formation (depth 0 or 1)
    from acyl chloride and amine.
    """
    late_amide_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal late_amide_formation_detected

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only consider late-stage reactions
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    reactants_part = rsmi.split(">")[0]
                    reactants = reactants_part.split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Checking reaction at depth {depth}: {rsmi}")

                    # Check if this is an acylation reaction
                    is_acylation = checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )

                    if is_acylation:
                        print(f"Found acylation reaction at depth {depth}")
                        late_amide_formation_detected = True
                        return

                    # If not detected by reaction type, check by functional groups
                    acyl_halide_present = False
                    amine_present = False
                    amide_in_product = False

                    # Check reactants for acyl halide and amine
                    for reactant in reactants:
                        if checker.check_fg("Acyl halide", reactant):
                            print(f"Found acyl halide in reactant: {reactant}")
                            acyl_halide_present = True

                        if checker.check_fg(
                            "Primary amine", reactant
                        ) or checker.check_fg("Secondary amine", reactant):
                            print(f"Found amine in reactant: {reactant}")
                            amine_present = True

                    # Check product for amide
                    if (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    ):
                        print(f"Found amide in product: {product}")
                        amide_in_product = True

                    # Check for Schotten-Baumann reaction specifically
                    is_schotten_baumann = checker.check_reaction(
                        "Schotten-Baumann to ester", rsmi
                    ) or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        rsmi,
                    )

                    if is_schotten_baumann:
                        print(f"Found Schotten-Baumann reaction at depth {depth}")
                        late_amide_formation_detected = True
                        return

                    # If we have all the components, it's likely an amide formation
                    if acyl_halide_present and amine_present and amide_in_product:
                        print(
                            f"Late-stage amide formation detected at depth {depth} by functional group analysis"
                        )
                        late_amide_formation_detected = True
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            if (
                not late_amide_formation_detected
            ):  # Stop traversal if we already found what we're looking for
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {late_amide_formation_detected}")
    return late_amide_formation_detected
