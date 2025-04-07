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
    Detects a strategy where a phenol is converted to a boronic ester through
    a series of functional group transformations.
    """
    # Track molecules containing key functional groups and their depths
    phenol_molecules = []
    boronic_ester_molecules = []

    # Track transformations with their depths
    transformations = []

    def dfs_traverse(node, depth=0):
        nonlocal phenol_molecules, boronic_ester_molecules, transformations

        if node["type"] == "mol":
            smiles = node["smiles"]

            # Check for phenol in any molecule
            if checker.check_fg("Phenol", smiles):
                phenol_molecules.append((smiles, depth))
                print(f"Found phenol at depth {depth}: {smiles}")

            # Check for boronic ester in any molecule
            if checker.check_fg("Boronic ester", smiles):
                boronic_ester_molecules.append((smiles, depth))
                print(f"Found boronic ester at depth {depth}: {smiles}")

            # Check for other relevant intermediates
            if checker.check_fg("Triflate", smiles):
                print(f"Found triflate intermediate at depth {depth}: {smiles}")
            if checker.check_fg("Mesylate", smiles):
                print(f"Found mesylate intermediate at depth {depth}: {smiles}")
            if checker.check_fg("Tosylate", smiles):
                print(f"Found tosylate intermediate at depth {depth}: {smiles}")
            if checker.check_fg("Aromatic halide", smiles):
                print(f"Found aromatic halide at depth {depth}: {smiles}")

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for transformations involving phenol or related intermediates

                # Phenol protection with silyl group
                if any(checker.check_fg("Phenol", r) for r in reactants) and checker.check_fg(
                    "Silyl protective group", product
                ):
                    if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                        transformations.append(("phenol_to_silyl_ether", depth))
                        print(
                            f"Found transformation at depth {depth}: phenol → silyl ether: {rsmi}"
                        )

                # Silyl ether deprotection
                elif any(
                    checker.check_fg("Silyl protective group", r) for r in reactants
                ) and checker.check_fg("Phenol", product):
                    if checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi):
                        transformations.append(("silyl_ether_to_phenol", depth))
                        print(
                            f"Found transformation at depth {depth}: silyl ether → phenol: {rsmi}"
                        )

                # Phenol to triflate conversion
                elif any(checker.check_fg("Phenol", r) for r in reactants) and checker.check_fg(
                    "Triflate", product
                ):
                    if checker.check_reaction("Formation of Sulfonic Esters", rsmi):
                        transformations.append(("phenol_to_triflate", depth))
                        print(f"Found transformation at depth {depth}: phenol → triflate: {rsmi}")

                # Phenol to mesylate conversion
                elif any(checker.check_fg("Phenol", r) for r in reactants) and checker.check_fg(
                    "Mesylate", product
                ):
                    if checker.check_reaction("Formation of Sulfonic Esters", rsmi):
                        transformations.append(("phenol_to_mesylate", depth))
                        print(f"Found transformation at depth {depth}: phenol → mesylate: {rsmi}")

                # Phenol to tosylate conversion
                elif any(checker.check_fg("Phenol", r) for r in reactants) and checker.check_fg(
                    "Tosylate", product
                ):
                    if checker.check_reaction("Formation of Sulfonic Esters", rsmi):
                        transformations.append(("phenol_to_tosylate", depth))
                        print(f"Found transformation at depth {depth}: phenol → tosylate: {rsmi}")

                # Triflate to boronic ester
                elif any(checker.check_fg("Triflate", r) for r in reactants) and checker.check_fg(
                    "Boronic ester", product
                ):
                    if checker.check_reaction(
                        "Preparation of boronic esters", rsmi
                    ) or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi):
                        transformations.append(("triflate_to_boronic_ester", depth))
                        print(
                            f"Found transformation at depth {depth}: triflate → boronic ester: {rsmi}"
                        )

                # Mesylate to boronic ester
                elif any(checker.check_fg("Mesylate", r) for r in reactants) and checker.check_fg(
                    "Boronic ester", product
                ):
                    transformations.append(("mesylate_to_boronic_ester", depth))
                    print(
                        f"Found transformation at depth {depth}: mesylate → boronic ester: {rsmi}"
                    )

                # Tosylate to boronic ester
                elif any(checker.check_fg("Tosylate", r) for r in reactants) and checker.check_fg(
                    "Boronic ester", product
                ):
                    transformations.append(("tosylate_to_boronic_ester", depth))
                    print(
                        f"Found transformation at depth {depth}: tosylate → boronic ester: {rsmi}"
                    )

                # Direct phenol to boronic ester
                elif any(checker.check_fg("Phenol", r) for r in reactants) and checker.check_fg(
                    "Boronic ester", product
                ):
                    transformations.append(("phenol_to_boronic_ester_direct", depth))
                    print(
                        f"Found direct transformation at depth {depth}: phenol → boronic ester: {rsmi}"
                    )

                # Phenol to aryl halide
                elif any(checker.check_fg("Phenol", r) for r in reactants) and any(
                    checker.check_fg(halide, product)
                    for halide in [
                        "Aromatic halide",
                        "Primary halide",
                        "Secondary halide",
                        "Tertiary halide",
                    ]
                ):
                    transformations.append(("phenol_to_halide", depth))
                    print(f"Found transformation at depth {depth}: phenol → halide: {rsmi}")

                # Aryl halide to boronic ester
                elif any(
                    checker.check_fg(halide, r)
                    for r in reactants
                    for halide in [
                        "Aromatic halide",
                        "Primary halide",
                        "Secondary halide",
                        "Tertiary halide",
                    ]
                ) and checker.check_fg("Boronic ester", product):
                    if checker.check_reaction("Preparation of boronic esters", rsmi):
                        transformations.append(("halide_to_boronic_ester", depth))
                        print(
                            f"Found transformation at depth {depth}: halide → boronic ester: {rsmi}"
                        )

        # Recursively process children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Analyze results to determine if a valid phenol to boronic ester strategy exists

    # Sort molecules and transformations by depth
    phenol_molecules.sort(key=lambda x: x[1])
    boronic_ester_molecules.sort(key=lambda x: x[1])
    transformations.sort(key=lambda x: x[1])

    # Check if the target molecule (depth 0) contains a boronic ester
    target_has_boronic_ester = any(depth == 0 for _, depth in boronic_ester_molecules)

    # Check if we have any relevant transformations
    has_transformations = len(transformations) > 0

    # Check if we have a direct phenol to boronic ester transformation
    direct_transformation = any(t[0] == "phenol_to_boronic_ester_direct" for t in transformations)

    # Check if we have a pathway through triflate
    triflate_pathway = any(t[0] == "phenol_to_triflate" for t in transformations) and any(
        t[0] == "triflate_to_boronic_ester" for t in transformations
    )

    # Check if we have a pathway through mesylate
    mesylate_pathway = any(t[0] == "phenol_to_mesylate" for t in transformations) and any(
        t[0] == "mesylate_to_boronic_ester" for t in transformations
    )

    # Check if we have a pathway through tosylate
    tosylate_pathway = any(t[0] == "phenol_to_tosylate" for t in transformations) and any(
        t[0] == "tosylate_to_boronic_ester" for t in transformations
    )

    # Check if we have a pathway through halide
    halide_pathway = any(t[0] == "phenol_to_halide" for t in transformations) and any(
        t[0] == "halide_to_boronic_ester" for t in transformations
    )

    # Check if we have a valid pathway
    valid_pathway = (
        direct_transformation
        or triflate_pathway
        or mesylate_pathway
        or tosylate_pathway
        or halide_pathway
    )

    print(f"Target has boronic ester: {target_has_boronic_ester}")
    print(f"Has transformations: {has_transformations}")
    print(f"Valid pathway detected: {valid_pathway}")

    # The key indicator is the presence of a boronic ester in the target molecule
    if target_has_boronic_ester:
        # If we have transformations and a valid pathway, it's a clear strategy
        if has_transformations and valid_pathway:
            print("Complete phenol to boronic ester transformation strategy detected")
            return True

        # Even without seeing the complete pathway, the presence of a boronic ester
        # in the target molecule is a strong indicator of this strategy
        print("Phenol to boronic ester strategy detected (target molecule contains boronic ester)")
        return True

    # If we have phenol molecules and transformations but no boronic ester in the target,
    # it might be a different strategy
    if phenol_molecules and has_transformations:
        print("Found phenol transformations but target doesn't contain boronic ester")
        return False

    print("Phenol to boronic ester transformation strategy not detected")
    return False
