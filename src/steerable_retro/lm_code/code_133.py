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
    Detects if the synthesis preserves a complex core structure throughout
    the route while modifying peripheral functional groups.
    """
    # Initialize variables
    molecules = []
    reactions = []

    # Track depth of each molecule in the synthesis tree
    molecule_depths = {}
    current_depth = 0

    def dfs_traverse(node, depth):
        nonlocal current_depth

        if node["type"] == "mol":
            molecules.append((node["smiles"], depth))
            molecule_depths[node["smiles"]] = depth

            # Track maximum depth
            current_depth = max(current_depth, depth)

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            reactions.append(node["metadata"]["rsmi"])

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Traverse the route to collect all molecules with their depths
    dfs_traverse(route, 0)

    if not molecules:
        print("No molecules found in the route")
        return False

    # Sort molecules by depth (ascending order - target molecule first)
    molecules.sort(key=lambda x: x[1])
    target_molecule = molecules[0][0]

    print(f"Target molecule: {target_molecule}")

    # Define potential core structures to check
    core_rings = [
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "thiophene",
        "furan",
        "pyrrole",
        "imidazole",
        "oxazole",
        "thiazole",
        "benzene",
        "naphthalene",
        "indole",
        "quinoline",
        "isoquinoline",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "cyclopropane",
        "cyclobutane",
        "cyclopentane",
        "cyclohexane",
    ]

    core_fgs = [
        "Sulfone",
        "Sulfonamide",
        "Carboxylic acid",
        "Ester",
        "Amide",
        "Urea",
        "Thiourea",
        "Nitro group",
        "Nitrile",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Ether",
        "Ketone",
        "Alcohol",
    ]

    # Check which cores are present in the target molecule
    target_cores = []
    for ring in core_rings:
        if checker.check_ring(ring, target_molecule):
            target_cores.append(("ring", ring))

    for fg in core_fgs:
        if checker.check_fg(fg, target_molecule):
            target_cores.append(("fg", fg))

    if not target_cores:
        print("No recognized core structures found in target molecule")
        return False

    print(f"Core structures in target molecule: {[core[1] for core in target_cores]}")

    # Check if core preservation strategy is used
    # We'll consider it a core preservation if:
    # 1. At least one core structure is preserved in most molecules (>60%)
    # 2. Most modifications occur in late-stage synthesis (last 30% of steps)

    core_preservation_counts = {core: 0 for core in target_cores}

    for mol_smiles, depth in molecules:
        for core_type, core_name in target_cores:
            if core_type == "ring" and checker.check_ring(core_name, mol_smiles):
                core_preservation_counts[(core_type, core_name)] += 1
            elif core_type == "fg" and checker.check_fg(core_name, mol_smiles):
                core_preservation_counts[(core_type, core_name)] += 1

    # Calculate preservation percentages
    total_molecules = len(molecules)
    preservation_percentages = {
        core: count / total_molecules * 100 for core, count in core_preservation_counts.items()
    }

    # Find the most preserved core
    most_preserved_core = max(preservation_percentages.items(), key=lambda x: x[1])

    print(f"Most preserved core: {most_preserved_core[0][1]} ({most_preserved_core[1]:.1f}%)")

    # Adjust threshold based on route length
    preservation_threshold = 60  # Default threshold
    if total_molecules <= 5:  # For short routes
        preservation_threshold = 50  # More lenient for short routes

    # Check if we have a significant core preservation
    significant_preservation = most_preserved_core[1] >= preservation_threshold

    if not significant_preservation:
        print(f"No significant core preservation detected (threshold: {preservation_threshold}%)")
        return False

    # Now check if modifications are primarily in late-stage synthesis
    # Get molecules that don't have the most preserved core
    modified_molecules = [
        (smiles, depth)
        for smiles, depth in molecules
        if (
            most_preserved_core[0][0] == "ring"
            and not checker.check_ring(most_preserved_core[0][1], smiles)
        )
        or (
            most_preserved_core[0][0] == "fg"
            and not checker.check_fg(most_preserved_core[0][1], smiles)
        )
    ]

    if not modified_molecules:
        # If all molecules have the core, it's technically core preservation
        print("Core is preserved in all molecules")
        return True

    # Calculate the threshold for "late-stage" (70% of max depth)
    late_stage_threshold = 0.7 * current_depth

    # Count how many modified molecules are in early vs late stage
    late_stage_mods = sum(1 for _, depth in modified_molecules if depth <= late_stage_threshold)
    early_stage_mods = len(modified_molecules) - late_stage_mods

    print(
        f"Modified molecules: {len(modified_molecules)} (Late stage: {late_stage_mods}, Early stage: {early_stage_mods})"
    )
    print(f"Late-stage threshold depth: {late_stage_threshold:.1f}")

    # Check if modifications are primarily late-stage (more than half)
    late_stage_modifications = late_stage_mods >= len(modified_molecules) / 2

    if not late_stage_modifications:
        print("Modifications are not primarily in late-stage synthesis")
        return False

    print("Core preservation strategy detected")
    return True
