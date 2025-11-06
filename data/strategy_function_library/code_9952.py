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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

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
    Detects if the synthesis follows a strategy of heterocyclic ring construction
    followed by late-stage nitro reduction.
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

    has_pyrazole_formation = False
    has_late_nitro_reduction = False
    max_depth = 0
    nitro_reduction_depth = None
    pyrazole_formation_depth = None

    def extract_depth(depth_info):
        # Extract numeric depth from the depth info string
        if not depth_info:
            return None
        match = re.search(r"Depth: (\d+)", depth_info)
        if match:
            return int(match.group(1))
        return None

    def dfs_traverse(node, current_depth=0):
        nonlocal has_pyrazole_formation, has_late_nitro_reduction
        nonlocal max_depth, nitro_reduction_depth, pyrazole_formation_depth
        nonlocal findings_json

        max_depth = max(max_depth, current_depth)

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Try to get depth from metadata, otherwise use traversal depth
                depth_info = node.get("metadata", {}).get("ID", "")
                metadata_depth = extract_depth(depth_info)
                node_depth = metadata_depth if metadata_depth is not None else current_depth

                print(f"Checking reaction at depth {node_depth}: {rsmi}")

                # Check for pyrazole formation
                if checker.check_ring("pyrazole", product) and not any(
                    checker.check_ring("pyrazole", r) for r in reactants
                ):
                    print(f"Found pyrazole formation at depth {node_depth}")
                    has_pyrazole_formation = True
                    pyrazole_formation_depth = node_depth
                    if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyrazole")
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                # Check for nitro reduction
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Found nitro reduction at depth {node_depth}")
                    has_late_nitro_reduction = True
                    nitro_reduction_depth = node_depth
                    if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

                # Fallback check for nitro reduction: nitro group disappears, amine appears
                if not has_late_nitro_reduction:
                    has_nitro_reactants = any(checker.check_fg("Nitro group", r) for r in reactants)
                    has_nitro_product = checker.check_fg("Nitro group", product)
                    has_amine_product = checker.check_fg(
                        "Primary amine", product
                    ) or checker.check_fg("Aniline", product)

                    if has_nitro_reactants and not has_nitro_product and has_amine_product:
                        print(f"Found nitro reduction (fallback check) at depth {node_depth}")
                        has_late_nitro_reduction = True
                        nitro_reduction_depth = node_depth
                        if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"] and checker.check_fg("Primary amine", product):
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        if "Aniline" not in findings_json["atomic_checks"]["functional_groups"] and checker.check_fg("Aniline", product):
                            findings_json["atomic_checks"]["functional_groups"].append("Aniline")
                        if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            if node["type"] == "reaction":
                # From reaction to chemical, depth remains the same
                dfs_traverse(child, current_depth)
            else:
                # From chemical to reaction, depth increases
                dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Synthesis analysis results:")
    print(f"Has pyrazole formation: {has_pyrazole_formation}, depth: {pyrazole_formation_depth}")
    print(f"Has nitro reduction: {has_late_nitro_reduction}, depth: {nitro_reduction_depth}")
    print(f"Max depth encountered: {max_depth}")

    result = False
    # Verify the sequence: pyrazole formation should happen at an earlier stage (higher depth)
    # than nitro reduction (lower depth)
    if has_pyrazole_formation and has_late_nitro_reduction:
        if pyrazole_formation_depth is not None and nitro_reduction_depth is not None:
            print(
                f"Comparing depths - Pyrazole formation: {pyrazole_formation_depth}, Nitro reduction: {nitro_reduction_depth}"
            )
            if pyrazole_formation_depth > nitro_reduction_depth:
                result = True
                # Add structural constraints if the condition is met
                findings_json["structural_constraints"].append({
                    "type": "co-occurrence",
                    "details": {
                        "targets": [
                            "ring_formation",
                            "Reduction of nitro groups to amines"
                        ]
                    }
                })
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "event_before": "ring_formation",
                        "target_before": "pyrazole",
                        "event_after": "Reduction of nitro groups to amines"
                    }
                })
        else:
            print("Cannot compare depths - one or both depths are None")

    return result, findings_json
