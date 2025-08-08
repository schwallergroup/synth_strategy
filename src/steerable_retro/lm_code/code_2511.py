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
    Detects if the synthesis route uses late-stage sulfonamide formation.
    """
    sulfonamide_formation_depths = []

    def dfs_traverse(node, current_depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a sulfonamide synthesis reaction using known reaction types
                is_sulfonamide_reaction = checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                ) or checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                )

                # Verify sulfonamide is formed (not present in reactants but present in product)
                sulfonamide_in_reactants = any(
                    checker.check_fg("Sulfonamide", r) for r in reactants
                )
                sulfonamide_in_product = checker.check_fg("Sulfonamide", product)

                # Check for the presence of required reactants
                sulfonyl_chloride_present = any(
                    checker.check_fg("Sulfonyl halide", r) for r in reactants
                )
                primary_amine_present = any(checker.check_fg("Primary amine", r) for r in reactants)
                secondary_amine_present = any(
                    checker.check_fg("Secondary amine", r) for r in reactants
                )

                # Determine if this is a sulfonamide formation reaction
                if is_sulfonamide_reaction or (
                    not sulfonamide_in_reactants
                    and sulfonamide_in_product
                    and sulfonyl_chloride_present
                    and (primary_amine_present or secondary_amine_present)
                ):
                    print(f"Found sulfonamide formation at depth: {current_depth}")
                    sulfonamide_formation_depths.append(current_depth)

        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)

    # Find the minimum depth where sulfonamide formation occurs
    min_depth = min(sulfonamide_formation_depths) if sulfonamide_formation_depths else None

    # Check if sulfonamide formation occurs at depth 0 or 1 (late stage)
    is_late_stage = min_depth is not None and min_depth <= 1
    print(f"Late-stage sulfonamide formation detected: {is_late_stage} (depth: {min_depth})")
    return is_late_stage
