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
    This function detects if halogens (F, Br) are retained throughout
    the synthesis until the final step.
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

    # Strategy JSON for reference (used to pick structural constraints)
    strategy_json_ref = {
      "function_id": "code_9165",
      "filepath": "../data/merged_good_perf/code_9165.py",
      "description": "This function detects if halogens (F, Br) are retained throughout the synthesis until the final step.",
      "atomic_checks": {
        "named_reactions": [],
        "ring_systems": [],
        "functional_groups": [
          "Aromatic halide",
          "Primary halide",
          "Secondary halide",
          "Tertiary halide",
          "Alkenyl halide"
        ]
      },
      "structural_constraints": [
        {
          "type": "co-occurrence",
          "details": {
            "targets": [
              "Fluorine in final product",
              "Fluorine in a starting material",
              "Bromine in a starting material"
            ],
            "description": "The strategy requires the co-occurrence of three conditions: Fluorine must be present in the final product, Fluorine must be present in at least one starting material, and Bromine must be present in at least one starting material."
          }
        },
        {
          "type": "positional",
          "details": {
            "target": "Fluorine-containing halide",
            "position": "last_stage",
            "description": "The final product molecule must contain a fluorine atom as part of a halide functional group."
          }
        },
        {
          "type": "positional",
          "details": {
            "target": "Fluorine-containing halide",
            "position": "first_stage",
            "description": "At least one starting material (a leaf node in the synthesis tree) must contain a fluorine atom as part of a halide functional group."
          }
        },
        {
          "type": "positional",
          "details": {
            "target": "Bromine-containing halide",
            "position": "first_stage",
            "description": "At least one starting material (a leaf node in the synthesis tree) must contain a bromine atom as part of a halide functional group."
          }
        }
      ]
    }

    # Track complete paths with halogen retention
    paths_with_f_retention = []
    paths_with_br_retention = []

    # Track current path during traversal
    current_path = []

    def dfs_traverse(node, depth=0, path_has_f=False, path_has_br=False):
        nonlocal paths_with_f_retention, paths_with_br_retention, current_path, findings_json
        # Add current node to path
        current_path.append(node)

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for halogens in this molecule
            has_f_aromatic = checker.check_fg("Aromatic halide", mol_smiles) and "F" in mol_smiles
            has_f_primary = checker.check_fg("Primary halide", mol_smiles) and "F" in mol_smiles
            has_f_secondary = checker.check_fg("Secondary halide", mol_smiles) and "F" in mol_smiles
            has_f_tertiary = checker.check_fg("Tertiary halide", mol_smiles) and "F" in mol_smiles
            has_f_alkenyl = checker.check_fg("Alkenyl halide", mol_smiles) and "F" in mol_smiles

            has_f = has_f_aromatic or has_f_primary or has_f_secondary or has_f_tertiary or has_f_alkenyl

            if has_f_aromatic and "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
            if has_f_primary and "Primary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Primary halide")
            if has_f_secondary and "Secondary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Secondary halide")
            if has_f_tertiary and "Tertiary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Tertiary halide")
            if has_f_alkenyl and "Alkenyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Alkenyl halide")

            has_br_aromatic = checker.check_fg("Aromatic halide", mol_smiles) and "Br" in mol_smiles
            has_br_primary = checker.check_fg("Primary halide", mol_smiles) and "Br" in mol_smiles
            has_br_secondary = checker.check_fg("Secondary halide", mol_smiles) and "Br" in mol_smiles
            has_br_tertiary = checker.check_fg("Tertiary halide", mol_smiles) and "Br" in mol_smiles
            has_br_alkenyl = checker.check_fg("Alkenyl halide", mol_smiles) and "Br" in mol_smiles

            has_br = has_br_aromatic or has_br_primary or has_br_secondary or has_br_tertiary or has_br_alkenyl

            if has_br_aromatic and "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
            if has_br_primary and "Primary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Primary halide")
            if has_br_secondary and "Secondary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Secondary halide")
            if has_br_tertiary and "Tertiary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Tertiary halide")
            if has_br_alkenyl and "Alkenyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Alkenyl halide")

            print(f"Molecule at depth {depth}: {mol_smiles}")
            print(f"  Has F: {has_f}, Has Br: {has_br}")

            # Update path halogen status
            path_has_f = path_has_f or has_f
            path_has_br = path_has_br or has_br

            # If this is a starting material (leaf node), check if path has retained halogens
            if not node.get("children", []):
                if path_has_f:
                    paths_with_f_retention.append(list(current_path))
                    # Record structural constraint: Fluorine in a starting material
                    constraint_found = next((c for c in strategy_json_ref["structural_constraints"] if c.get("details", {}).get("target") == "Fluorine-containing halide" and c.get("details", {}).get("position") == "first_stage"), None)
                    if constraint_found and constraint_found not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint_found)
                    print(f"  Complete path with F retention found at depth {depth}")

                if path_has_br:
                    paths_with_br_retention.append(list(current_path))
                    # Record structural constraint: Bromine in a starting material
                    constraint_found = next((c for c in strategy_json_ref["structural_constraints"] if c.get("details", {}).get("target") == "Bromine-containing halide" and c.get("details", {}).get("position") == "first_stage"), None)
                    if constraint_found and constraint_found not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint_found)
                    print(f"  Complete path with Br retention found at depth {depth}")

        # Process children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # Depth increases only when going from chemical to reaction
                new_depth = depth + 1
            dfs_traverse(child, new_depth, path_has_f, path_has_br)

        # Remove current node from path when backtracking
        current_path.pop()

    # Start traversal
    dfs_traverse(route)

    # Check if the final product (depth 0) has F or Br
    final_product_smiles = route["smiles"]
    final_has_f_aromatic = checker.check_fg("Aromatic halide", final_product_smiles) and "F" in final_product_smiles
    final_has_f_primary = checker.check_fg("Primary halide", final_product_smiles) and "F" in final_product_smiles
    final_has_f_secondary = checker.check_fg("Secondary halide", final_product_smiles) and "F" in final_product_smiles
    final_has_f_tertiary = checker.check_fg("Tertiary halide", final_product_smiles) and "F" in final_product_smiles
    final_has_f_alkenyl = checker.check_fg("Alkenyl halide", final_product_smiles) and "F" in final_product_smiles

    final_has_f = final_has_f_aromatic or final_has_f_primary or final_has_f_secondary or final_has_f_tertiary or final_has_f_alkenyl

    if final_has_f_aromatic and "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
        findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
    if final_has_f_primary and "Primary halide" not in findings_json["atomic_checks"]["functional_groups"]:
        findings_json["atomic_checks"]["functional_groups"].append("Primary halide")
    if final_has_f_secondary and "Secondary halide" not in findings_json["atomic_checks"]["functional_groups"]:
        findings_json["atomic_checks"]["functional_groups"].append("Secondary halide")
    if final_has_f_tertiary and "Tertiary halide" not in findings_json["atomic_checks"]["functional_groups"]:
        findings_json["atomic_checks"]["functional_groups"].append("Tertiary halide")
    if final_has_f_alkenyl and "Alkenyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
        findings_json["atomic_checks"]["functional_groups"].append("Alkenyl halide")

    final_has_br_aromatic = checker.check_fg("Aromatic halide", final_product_smiles) and "Br" in final_product_smiles
    final_has_br_primary = checker.check_fg("Primary halide", final_product_smiles) and "Br" in final_product_smiles
    final_has_br_secondary = checker.check_fg("Secondary halide", final_product_smiles) and "Br" in final_product_smiles
    final_has_br_tertiary = checker.check_fg("Tertiary halide", final_product_smiles) and "Br" in final_product_smiles
    final_has_br_alkenyl = checker.check_fg("Alkenyl halide", final_product_smiles) and "Br" in final_product_smiles

    final_has_br = final_has_br_aromatic or final_has_br_primary or final_has_br_secondary or final_has_br_tertiary or final_has_br_alkenyl

    if final_has_br_aromatic and "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
        findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
    if final_has_br_primary and "Primary halide" not in findings_json["atomic_checks"]["functional_groups"]:
        findings_json["atomic_checks"]["functional_groups"].append("Primary halide")
    if final_has_br_secondary and "Secondary halide" not in findings_json["atomic_checks"]["functional_groups"]:
        findings_json["atomic_checks"]["functional_groups"].append("Secondary halide")
    if final_has_br_tertiary and "Tertiary halide" not in findings_json["atomic_checks"]["functional_groups"]:
        findings_json["atomic_checks"]["functional_groups"].append("Tertiary halide")
    if final_has_br_alkenyl and "Alkenyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
        findings_json["atomic_checks"]["functional_groups"].append("Alkenyl halide")

    print(f"Final product has F: {final_has_f}, has Br: {final_has_br}")
    print(f"Paths with F retention: {len(paths_with_f_retention)}")
    print(f"Paths with Br retention: {len(paths_with_br_retention)}")

    # Halogen retention strategy is detected if:
    # 1. Final product has F and there's at least one path with F retention
    # 2. There's at least one path with Br retention (Br doesn't need to be in final product)
    f_strategy = final_has_f and len(paths_with_f_retention) > 0
    br_strategy = len(paths_with_br_retention) > 0

    if final_has_f:
        # Record structural constraint: Fluorine in final product
        constraint_found = next((c for c in strategy_json_ref["structural_constraints"] if c.get("details", {}).get("target") == "Fluorine-containing halide" and c.get("details", {}).get("position") == "last_stage"), None)
        if constraint_found and constraint_found not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(constraint_found)

    strategy_found = f_strategy and br_strategy

    if strategy_found:
        # Record co-occurrence constraint if all conditions are met
        co_occurrence_constraint = next((c for c in strategy_json_ref["structural_constraints"] if c.get("type") == "co-occurrence"), None)
        if co_occurrence_constraint and co_occurrence_constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(co_occurrence_constraint)

    print(f"Halogen retention strategy detected: {strategy_found}")
    return strategy_found, findings_json