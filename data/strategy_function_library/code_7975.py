from typing import Tuple, Dict, List
import copy
from rdkit.Chem import AllChem, rdFMCS
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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a sequence of nitrogen oxidation state changes:
    nitro → amine → hydrazine → heterocycle (pyrazole)

    In retrosynthetic traversal, we expect to find these in reverse order:
    pyrazole (target) → hydrazine → amine → nitro (starting material)
    """
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    strategy_json = {
      "function_id": "code_7975",
      "filepath": "../data/merged_good_perf/code_7975.py",
      "description": "Detects a specific sequence of nitrogen oxidation state changes. In the forward synthesis direction, this is: nitro -> amine -> hydrazine -> pyrazole. The function validates this by checking the reverse sequence in a retrosynthetic route: pyrazole formation, preceded by hydrazine synthesis, preceded by amine formation from a nitro group.",
      "atomic_checks": {
        "named_reactions": [
          "Reduction of nitro groups to amines",
          "Hydrazine synthesis from amine",
          "pyrazole",
          "Pyrazole formation",
          "ring_formation"
        ],
        "ring_systems": [
          "pyrazole"
        ],
        "functional_groups": [
          "Nitro group",
          "Primary amine",
          "Hydrazine"
        ]
      },
      "structural_constraints": [
        {
          "type": "co-occurrence",
          "details": {
            "targets": [
              "Nitro group",
              "Reduction of nitro groups to amines",
              "Hydrazine synthesis from amine",
              "Pyrazole formation"
            ],
            "description": "The route must contain all four key chemical entities/transformations, tracked by the flags: nitro_found, amine_found, hydrazine_found, and pyrazole_found."
          }
        },
        {
          "type": "sequence",
          "details": {
            "ordered_events": [
              "Pyrazole formation",
              "Hydrazine synthesis from amine",
              "Reduction of nitro groups to amines",
              "Nitro group"
            ],
            "relation": "depth_le",
            "description": "The events must occur in a specific order in the retrosynthetic tree, checked via `pyrazole_depth <= hydrazine_depth <= amine_depth <= nitro_depth`. The pyrazole formation must be closest to the target (lowest depth), and the nitro group must be closest to the starting materials (highest depth)."
          }
        }
      ]
    }

    # Track the stages we've found
    nitro_found = False
    amine_found = False
    hydrazine_found = False
    pyrazole_found = False

    # Track the order (depth) of transformations
    # In retrosynthesis, lower depth = later stage (closer to target)
    nitro_depth = -1
    amine_depth = -1
    hydrazine_depth = -1
    pyrazole_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal nitro_found, amine_found, hydrazine_found, pyrazole_found
        nonlocal nitro_depth, amine_depth, hydrazine_depth, pyrazole_depth
        nonlocal findings_json

        if node["type"] == "mol":
            # Check for functional groups in molecules
            mol_smiles = node["smiles"]

            # Record the occurrence of each group, updating if found at an earlier stage
            if checker.check_fg("Nitro group", mol_smiles) and (
                nitro_depth == -1 or depth < nitro_depth
            ):
                print(f"Found nitro group in molecule at depth {depth}")
                nitro_found = True
                nitro_depth = depth
                if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

        elif node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro reduction to amine
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Found nitro reduction to amine at depth {depth}")
                    if amine_depth == -1 or depth < amine_depth:
                        amine_found = True
                        amine_depth = depth
                        if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

                # Check for amine to hydrazine transformation
                if checker.check_reaction("Hydrazine synthesis from amine", rsmi):
                    print(
                        f"Found amine to hydrazine transformation (specific reaction) at depth {depth}"
                    )
                    if hydrazine_depth == -1 or depth < hydrazine_depth:
                        hydrazine_found = True
                        hydrazine_depth = depth
                        if "Hydrazine synthesis from amine" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Hydrazine synthesis from amine")

                # Generic check for amine to hydrazine transformation
                if (
                    any(checker.check_fg("Primary amine", r) for r in reactants)
                    and checker.check_fg("Hydrazine", product)
                    and not any(checker.check_fg("Hydrazine", r) for r in reactants)
                ):
                    print(f"Found amine to hydrazine transformation at depth {depth}")
                    if hydrazine_depth == -1 or depth < hydrazine_depth:
                        hydrazine_found = True
                        hydrazine_depth = depth
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        if "Hydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Hydrazine")

                # Check for hydrazine to pyrazole transformation
                if checker.check_reaction("pyrazole", rsmi) or checker.check_reaction(
                    "Pyrazole formation", rsmi
                ):
                    print(f"Found pyrazole formation reaction at depth {depth}")
                    if pyrazole_depth == -1 or depth < pyrazole_depth:
                        pyrazole_found = True
                        pyrazole_depth = depth
                        if "pyrazole" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("pyrazole")
                        if "Pyrazole formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Pyrazole formation")

                # Generic check for hydrazine to pyrazole transformation
                if (
                    any(checker.check_fg("Hydrazine", r) for r in reactants)
                    and checker.check_ring("pyrazole", product)
                    and not any(checker.check_ring("pyrazole", r) for r in reactants)
                ):
                    print(f"Found hydrazine to pyrazole transformation at depth {depth}")
                    if pyrazole_depth == -1 or depth < pyrazole_depth:
                        pyrazole_found = True
                        pyrazole_depth = depth
                        if "Hydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Hydrazine")
                        if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("pyrazole")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth when going from chemical to reaction
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitro found: {nitro_found} at depth {nitro_depth}")
    print(f"Amine found: {amine_found} at depth {amine_depth}")
    print(f"Hydrazine found: {hydrazine_found} at depth {hydrazine_depth}")
    print(f"Pyrazole found: {pyrazole_found} at depth {pyrazole_depth}")

    # Check if we found all stages in the correct order for retrosynthesis
    # In retrosynthetic traversal, target (pyrazole) has lowest depth, starting material (nitro) has highest
    # Allow for equal depths in case multiple functional groups are in the same molecule
    correct_sequence = (
        nitro_found
        and amine_found
        and hydrazine_found
        and pyrazole_found
        and pyrazole_depth <= hydrazine_depth <= amine_depth <= nitro_depth
    )

    if nitro_found and amine_found and hydrazine_found and pyrazole_found:
        # Add the co-occurrence constraint if all flags are true
        co_occurrence_constraint = next((c for c in strategy_json["structural_constraints"] if c["type"] == "co-occurrence"), None)
        if co_occurrence_constraint and co_occurrence_constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(co_occurrence_constraint)

    if correct_sequence:
        # Add the sequence constraint if the order is correct
        sequence_constraint = next((c for c in strategy_json["structural_constraints"] if c["type"] == "sequence"), None)
        if sequence_constraint and sequence_constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(sequence_constraint)

    print(f"Correct sequence: {correct_sequence}")
    return correct_sequence, findings_json