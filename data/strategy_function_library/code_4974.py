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


BOC_PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of primary amine",
    "Boc amine protection of secondary amine",
]

BOC_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a strategy involving early-stage Boc amine protection. Early-stage is defined as occurring in the first half of the synthesis sequence (i.e., at a depth greater than the midpoint of the total synthesis depth). The function identifies specific Boc protection reactions from the BOC_PROTECTION_REACTIONS list and ignores routes where the Boc group is present on a starting material.
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

    # Track protection and deprotection events with their depths
    protection_events = []
    deprotection_events = []
    min_depth = float("inf")
    max_depth = 0
    result = False

    # Get the target molecule (root of the tree)
    target_mol = route["smiles"]
    print(f"Target molecule: {target_mol}")

    # Check if target molecule contains a Boc-protected amine
    target_has_boc = checker.check_fg("Boc", target_mol)
    if target_has_boc:
        print(f"Target molecule contains Boc group")
        if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
            findings_json["atomic_checks"]["functional_groups"].append("Boc")

    def dfs_traverse(node, current_depth=0):
        nonlocal min_depth, max_depth, findings_json

        # Update depth tracking for all nodes
        if node["type"] == "reaction":
            min_depth = min(min_depth, current_depth)
            max_depth = max(max_depth, current_depth)

            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for amine protection reactions
                is_protection = False
                for r in BOC_PROTECTION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        is_protection = True
                        protection_events.append((current_depth, rsmi))
                        print(f"Found Boc protection at depth: {current_depth}")
                        if r not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)
                        break

                # Check for deprotection reactions
                is_deprotection = False
                for r in BOC_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        is_deprotection = True
                        deprotection_events.append((current_depth, rsmi))
                        print(f"Found Boc deprotection at depth: {current_depth}")
                        if r not in findings_json["atomic_checks"]["named_reactions"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)
                        break
        elif node["type"] == "mol":
            # Track depth for molecule nodes too
            min_depth = min(min_depth, current_depth)
            max_depth = max(max_depth, current_depth)

            # Check if this is a starting material with Boc
            if node.get("in_stock", False) and checker.check_fg("Boc", node["smiles"]):
                print(f"Found starting material with Boc at depth: {current_depth}")
                if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boc")

        # Continue traversal
        for child in node.get("children", []):
            # New logic: depth increases only from chemical to reaction node
            if node["type"] == "reaction":
                dfs_traverse(child, current_depth)
            else: # node["type"] == "mol"
                dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # If no protection found and target doesn't have Boc, or only one depth level, return False
    if (not protection_events and not target_has_boc) or max_depth == min_depth:
        print("No amine protection found or synthesis has only one depth level")
        result = False
        return result, findings_json

    # Calculate the midpoint of the synthesis (in retrosynthetic direction)
    midpoint = min_depth + (max_depth - min_depth) / 2
    print(f"Synthesis depth range: {min_depth} to {max_depth}, Midpoint: {midpoint}")

    # If target has Boc but no protection events, it might be from starting materials
    if target_has_boc and not protection_events:
        print("Target has Boc group but no protection reaction found")
        # Check if any starting material has Boc
        has_boc_starting_material = False

        def check_starting_materials(node):
            nonlocal has_boc_starting_material, findings_json
            if node["type"] == "mol" and node.get("in_stock", False):
                if checker.check_fg("Boc", node["smiles"]):
                    has_boc_starting_material = True
                    # Record the structural constraint if a starting material has Boc
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "Boc",
                            "scope": "starting_material",
                            "description": "The route is invalid if any starting material contains a Boc functional group."
                        }
                    })
            for child in node.get("children", []):
                check_starting_materials(child)

        check_starting_materials(route)

        if has_boc_starting_material:
            print("Boc group comes from starting materials, not a protection strategy")
            result = False
            return result, findings_json
        else:
            print("Boc group must be introduced in synthesis but protection reaction not detected")
            result = True
            return result, findings_json

    # In retrosynthesis, late-stage protection means protection depth is GREATER than midpoint
    # (protection happens early in retrosynthesis = late in forward synthesis)
    late_stage_protections = [depth for depth, _ in protection_events if depth >= midpoint]

    if late_stage_protections:
        print(f"Late-stage protection depths: {late_stage_protections}")
        print("Late-stage amine protection strategy detected")
        result = True
        # Record the structural constraint for early-stage protection (which is late-stage in retrosynthesis)
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Boc amine protection",
                "position": "first_half",
                "description": "A Boc amine protection reaction must occur in the first half of the forward synthesis route (retrosynthetic depth >= midpoint)."
            }
        })
        return result, findings_json

    print(f"Protection depths: {[d for d, _ in protection_events]}")
    print("No late-stage amine protection strategy detected")
    result = False
    return result, findings_json
