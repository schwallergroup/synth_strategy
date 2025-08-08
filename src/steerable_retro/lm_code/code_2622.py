#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
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
    This function detects if the synthetic route employs a late-stage nitro reduction
    as the final step to generate an amine functional group.
    """
    # Track if we found a nitro reduction in the final step
    found_late_stage_nitro_reduction = False

    print("Starting late_stage_nitro_reduction_strategy analysis...")

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_nitro_reduction

        print(f"Examining node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and depth <= 1:  # Final or near-final step
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing reaction SMILES at depth {depth}: {rsmi}")

                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains amine group (any type)
                has_amine = (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                    or checker.check_fg("Aniline", product)
                )

                if has_amine:
                    print(f"Product contains amine: {product}")

                    # Check if reactant contains nitro group and the reaction is a nitro reduction
                    for reactant in reactants:
                        if checker.check_fg("Nitro group", reactant):
                            print(f"Found reactant with nitro group: {reactant}")

                            # First try the specific reaction check
                            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                                print(
                                    f"CONFIRMED: Found nitro reduction in late stage (depth {depth}): {rsmi}"
                                )
                                found_late_stage_nitro_reduction = True
                                break
                            # If that fails, check if this is a general reduction with nitro -> amine conversion
                            else:
                                print(f"Checking atom mapping for nitro reduction...")
                                # Get atom indices to verify the nitro group is converted to amine
                                try:
                                    nitro_indices = checker.get_fg_atom_indices(
                                        "Nitro group", reactant
                                    )
                                    if nitro_indices and len(nitro_indices) > 0:
                                        # The first atom in nitro group is typically the N atom
                                        for nitro_group in nitro_indices:
                                            for atom_tuple in nitro_group:
                                                # Look for the nitrogen atom in the nitro group
                                                for atom_idx in atom_tuple:
                                                    # Check if this atom is mapped to an amine in the product
                                                    if (
                                                        f"[NH2:{atom_idx}]" in product
                                                        or f"NH2:{atom_idx}" in product
                                                    ):
                                                        print(
                                                            f"CONFIRMED: Found nitro reduction based on atom mapping (depth {depth}): {rsmi}"
                                                        )
                                                        found_late_stage_nitro_reduction = True
                                                        break
                                                if found_late_stage_nitro_reduction:
                                                    break
                                            if found_late_stage_nitro_reduction:
                                                break
                                except Exception as e:
                                    print(f"Error checking atom mapping: {e}")

                                # If we still haven't found it, look for common reagents for nitro reduction
                                if not found_late_stage_nitro_reduction:
                                    reagents = rsmi.split(">")[1].split(".")
                                    if any(r for r in reagents if "[Fe]" in r or "Fe" in r):
                                        print(
                                            f"Found iron reagent, likely a nitro reduction: {rsmi}"
                                        )
                                        found_late_stage_nitro_reduction = True
                                        break
                        else:
                            print(f"Reactant does not contain nitro group: {reactant}")
                else:
                    print(f"Product does not contain any amine group: {product}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    print(f"Final result: found_late_stage_nitro_reduction = {found_late_stage_nitro_reduction}")
    return found_late_stage_nitro_reduction
