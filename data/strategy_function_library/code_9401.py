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
    Detects a late-stage urea or thiourea formation strategy. It identifies specific
    named reactions (e.g., from isocyanates) and also general cases where a urea or
    thiourea functional group is formed in one of the final synthetic steps.
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

    # Initialize tracking variables
    has_urea_formation = False
    reaction_depths = {}
    max_depth = 0
    urea_formation_depths = []
    amide_formation_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal has_urea_formation, max_depth, findings_json

        # Track maximum depth to determine early vs late stage
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Store reaction at its depth for later analysis
            if depth not in reaction_depths:
                reaction_depths[depth] = []
            reaction_depths[depth].append(node)

            # Check if this reaction is a urea formation
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for urea formation reactions
                urea_reactions = [
                    "Urea synthesis via isocyanate and primary amine",
                    "Urea synthesis via isocyanate and secondary amine",
                    "Urea synthesis via isocyanate and diazo",
                    "Urea synthesis via isocyanate and sulfonamide"
                ]
                detected_named_reaction = False
                for r_name in urea_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        if (checker.check_fg("Urea", product) or checker.check_fg("Thiourea", product)):
                            print(f"Detected urea formation at depth {depth}")
                            has_urea_formation = True
                            urea_formation_depths.append(depth)
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                            if checker.check_fg("Urea", product):
                                findings_json["atomic_checks"]["functional_groups"].append("Urea")
                            if checker.check_fg("Thiourea", product):
                                findings_json["atomic_checks"]["functional_groups"].append("Thiourea")
                            detected_named_reaction = True
                            break

                # Also check for urea formation by looking at functional groups
                if not detected_named_reaction and \
                   (checker.check_fg("Urea", product) or checker.check_fg("Thiourea", product)) and \
                   not any(
                       checker.check_fg("Urea", r) or checker.check_fg("Thiourea", r)
                       for r in reactants
                   ):
                    print(f"Detected urea formation (by FG analysis) at depth {depth}")
                    has_urea_formation = True
                    urea_formation_depths.append(depth)
                    if checker.check_fg("Urea", product):
                        findings_json["atomic_checks"]["functional_groups"].append("Urea")
                        findings_json["atomic_checks"]["named_reactions"].append("Urea_formation")
                    if checker.check_fg("Thiourea", product):
                        findings_json["atomic_checks"]["functional_groups"].append("Thiourea")
                        findings_json["atomic_checks"]["named_reactions"].append("Thiourea_formation")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Urea formation depths: {urea_formation_depths}")
    print(f"Amide formation depths: {amide_formation_depths}")

    # Define late stage as the first third of the synthesis
    late_stage_threshold = max(1, max_depth // 3)
    print(f"Max depth: {max_depth}, Late stage threshold: {late_stage_threshold}")

    # Check if urea formation occurs in the late stage
    late_stage_urea = any(depth <= late_stage_threshold for depth in urea_formation_depths)
    if late_stage_urea:
        print(f"Urea formation occurs at late stage (depth <= {late_stage_threshold})")
        # Add the structural constraint if late stage urea formation is detected
        # This assumes the 'positional' constraint is met if any of the target reactions/FGs are found late stage
        if has_urea_formation: # Ensure a urea formation was actually found
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "targets": [
                        "Urea synthesis via isocyanate and primary amine",
                        "Urea synthesis via isocyanate and secondary amine",
                        "Urea synthesis via isocyanate and diazo",
                        "Urea synthesis via isocyanate and sulfonamide",
                        "Urea_formation",
                        "Thiourea_formation"
                    ],
                    "position": "late_stage"
                }
            })

    # Check if there's an amide formation/disconnection in an earlier stage
    early_stage_amide = any(depth > late_stage_threshold for depth in amide_formation_depths)
    if early_stage_amide:
        print(
            f"Amide formation/disconnection occurs at early stage (depth > {late_stage_threshold})"
        )

    # The strategy is present if we have a late-stage urea formation and an earlier amide formation/disconnection
    result = has_urea_formation and (late_stage_urea or early_stage_amide)
    print(f"Late-stage urea formation strategy detected: {result}")
    return result, findings_json
