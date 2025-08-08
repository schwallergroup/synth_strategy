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
    Detects if the synthesis route uses a triflate as a leaving group
    for nucleophilic substitution.
    """
    # Track reaction nodes with depth information
    triflate_formation_reactions = []
    triflate_displacement_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for triflate formation
                if checker.check_reaction("Alcohol to triflate conversion", rsmi):
                    print(f"Detected triflate formation reaction at depth {depth}: {rsmi}")
                    triflate_formation_reactions.append((depth, rsmi))
                elif any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    or checker.check_fg("Phenol", r)
                    for r in reactants
                ) and checker.check_fg("Triflate", product):
                    print(f"Detected triflate formation at depth {depth}: {rsmi}")
                    triflate_formation_reactions.append((depth, rsmi))

                # Check for triflate displacement
                if any(checker.check_fg("Triflate", r) for r in reactants) and not checker.check_fg(
                    "Triflate", product
                ):
                    # Check for common nucleophilic substitution reactions
                    if (
                        checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                        or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                        or checker.check_reaction("Sonogashira alkyne_aryl OTf", rsmi)
                        or checker.check_reaction("Sonogashira acetylene_aryl OTf", rsmi)
                        or checker.check_reaction("Sonogashira alkyne_alkenyl OTf", rsmi)
                        or checker.check_reaction("Sonogashira acetylene_alkenyl OTf", rsmi)
                        or checker.check_reaction("Stille reaction_vinyl OTf", rsmi)
                        or checker.check_reaction("Stille reaction_aryl OTf", rsmi)
                        or checker.check_reaction("Stille reaction_benzyl OTf", rsmi)
                        or checker.check_reaction("Stille reaction_allyl OTf", rsmi)
                        or checker.check_reaction("Stille reaction_other OTf", rsmi)
                        or checker.check_reaction(
                            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                        )
                    ):
                        print(
                            f"Detected triflate displacement via coupling reaction at depth {depth}: {rsmi}"
                        )
                        triflate_displacement_reactions.append((depth, rsmi))
                    # Check for nucleophilic substitution with amines
                    elif any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Tertiary amine", r)
                        for r in reactants
                    ) and (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                    ):
                        print(f"Detected triflate displacement by amine at depth {depth}: {rsmi}")
                        triflate_displacement_reactions.append((depth, rsmi))
                    # Check for other nucleophilic substitutions
                    elif any(
                        checker.check_fg("Phenol", r)
                        or checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        for r in reactants
                    ):
                        print(f"Detected triflate displacement by alcohol at depth {depth}: {rsmi}")
                        triflate_displacement_reactions.append((depth, rsmi))

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have both formation and displacement
    has_triflate_formation = len(triflate_formation_reactions) > 0
    has_triflate_displacement = len(triflate_displacement_reactions) > 0

    # Check if triflate is used as a starting material
    def check_starting_triflate():
        def check_mol_node(node):
            if node["type"] == "mol" and node.get("in_stock", False):
                if checker.check_fg("Triflate", node["smiles"]):
                    return True
            for child in node.get("children", []):
                if check_mol_node(child):
                    return True
            return False

        return check_mol_node(route)

    has_starting_triflate = check_starting_triflate()

    # In retrosynthetic direction, formation should have higher depth than displacement
    # In forward synthesis, we form triflate first, then displace it
    def check_sequence():
        if not triflate_formation_reactions or not triflate_displacement_reactions:
            return False

        # Get minimum depth for each type of reaction
        # Lower depth = later stage in synthesis (closer to target)
        min_formation_depth = min([depth for depth, _ in triflate_formation_reactions])
        min_displacement_depth = min([depth for depth, _ in triflate_displacement_reactions])

        # In retrosynthesis, formation should be at a lower depth (later stage)
        # than displacement (earlier stage)
        print(
            f"Min formation depth: {min_formation_depth}, Min displacement depth: {min_displacement_depth}"
        )
        return min_formation_depth <= min_displacement_depth

    sequence_correct = check_sequence()

    # If we have a starting triflate, we don't need formation
    result = (has_triflate_formation or has_starting_triflate) and has_triflate_displacement

    # If we have both formation and displacement, check sequence
    if has_triflate_formation and has_triflate_displacement:
        result = result and sequence_correct

    print(
        f"Triflate formation: {has_triflate_formation}, Triflate displacement: {has_triflate_displacement}"
    )
    print(f"Starting triflate: {has_starting_triflate}, Sequence correct: {sequence_correct}")
    print(f"Final result: {result}")
    return result
