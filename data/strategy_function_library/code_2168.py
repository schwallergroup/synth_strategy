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
    This function detects a strategy where a silyl protection group (TBDMS) is installed
    early and maintained throughout most of the synthesis.
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
    strategy_json = {
      "function_id": "code_2168",
      "filepath": "../data/merged_good_perf/code_2168.py",
      "description": "This function detects a strategy where a silyl protection group (TBDMS) is installed and maintained throughout most of the synthesis.",
      "atomic_checks": {
        "named_reactions": [
          "Alcohol protection with silyl ethers"
        ],
        "ring_systems": [],
        "functional_groups": [
          "Silyl protective group"
        ]
      },
      "structural_constraints": [
        {
          "type": "positional",
          "details": {
            "target": "Silyl protective group",
            "position": "late_stage"
          }
        },
        {
          "type": "co-occurrence",
          "details": {
            "targets": [
              "Alcohol protection with silyl ethers",
              "Silyl protective group"
            ]
          }
        },
        {
          "type": "sequence",
          "details": {
            "target": "Silyl protective group",
            "start_position": "starting_material",
            "end_position": "late_stage"
          }
        },
        {
          "type": "count",
          "details": {
            "target": "synthesis_stages_with_silyl_group",
            "operator": ">=",
            "value": 2
          }
        }
      ]
    }

    protection_step_found = False
    protection_maintained = False
    protected_mol_at_depth = {}  # Track protected molecules at each depth

    # Check if starting material is already protected
    def check_starting_materials(node):
        if node["type"] == "mol" and node.get("in_stock", False):
            if checker.check_fg("Silyl protective group", node["smiles"]):
                findings_json["atomic_checks"]["functional_groups"].append("Silyl protective group")
                return True
        for child in node.get("children", []):
            if check_starting_materials(child):
                return True
        return False

    starting_material_protected = check_starting_materials(route)
    if starting_material_protected:
        print("Starting material already contains silyl protection")

    def dfs_traverse(node, current_depth=0):
        nonlocal protection_step_found, protection_maintained

        node_depth = current_depth

        if node["type"] == "mol":
            # Check if molecule has silyl protection
            mol_smiles = node["smiles"]
            if checker.check_fg("Silyl protective group", mol_smiles):
                if "Silyl protective group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Silyl protective group")
                protected_mol_at_depth[node_depth] = mol_smiles
                print(f"Detected protected molecule at depth {node_depth}")

                # If we see protection at a late stage (low depth)
                if node_depth <= 3:
                    protection_maintained = True
                    # Add positional constraint if met
                    for constraint in strategy_json["structural_constraints"]:
                        if constraint["type"] == "positional" and constraint["details"]["target"] == "Silyl protective group" and constraint["details"]["position"] == "late_stage":
                            if constraint not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(constraint)
                    print(f"Confirmed silyl protection at late stage (depth {node_depth})")

        elif node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product contains silyl protected oxygen
            product_has_silyl = checker.check_fg("Silyl protective group", product_smiles)
            if product_has_silyl and "Silyl protective group" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Silyl protective group")

            # Check if any reactant already has silyl protection
            reactants_have_silyl = any(
                checker.check_fg("Silyl protective group", r)
                for r in reactants_smiles
            )
            if reactants_have_silyl and "Silyl protective group" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Silyl protective group")

            # Check if this is a silyl protection reaction
            is_protection_reaction = checker.check_reaction(
                "Alcohol protection with silyl ethers", rsmi
            )
            if is_protection_reaction:
                if "Alcohol protection with silyl ethers" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Alcohol protection with silyl ethers")

            # If this is a protection step (early in synthesis)
            if is_protection_reaction or (product_has_silyl and not reactants_have_silyl):
                protection_step_found = True
                protected_mol_at_depth[node_depth] = product_smiles
                # Add co-occurrence constraint if met
                for constraint in strategy_json["structural_constraints"]:
                    if constraint["type"] == "co-occurrence" and "Alcohol protection with silyl ethers" in constraint["details"]["targets"] and "Silyl protective group" in constraint["details"]["targets"]:
                        if constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(constraint)
                print(f"Detected silyl protection installation at depth {node_depth}")

            # If protection is maintained in subsequent steps
            elif product_has_silyl and reactants_have_silyl:
                protected_mol_at_depth[node_depth] = product_smiles
                print(f"Detected silyl protection maintained at depth {node_depth}")

        # Continue traversing
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, node_depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, node_depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Verify we have evidence of protection strategy
    has_protection_strategy = False

    # Case 1: Protection step found and maintained to late stages
    if protection_step_found and protection_maintained:
        depths = sorted(protected_mol_at_depth.keys())
        if len(depths) >= 2:
            has_protection_strategy = True
            # Add sequence constraint if met
            for constraint in strategy_json["structural_constraints"]:
                if constraint["type"] == "sequence" and constraint["details"]["target"] == "Silyl protective group" and constraint["details"]["start_position"] == "starting_material" and constraint["details"]["end_position"] == "late_stage":
                    if constraint not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint)
            # Add count constraint if met
            for constraint in strategy_json["structural_constraints"]:
                if constraint["type"] == "count" and constraint["details"]["target"] == "synthesis_stages_with_silyl_group" and constraint["details"]["operator"] == ">=" and constraint["details"]["value"] == 2:
                    if constraint not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint)
            print(f"Detected silyl protection strategy from depth {max(depths)} to {min(depths)}")

    # Case 2: Starting material was already protected and maintained
    elif starting_material_protected and protection_maintained:
        depths = sorted(protected_mol_at_depth.keys())
        if len(depths) >= 2:
            has_protection_strategy = True
            # Add sequence constraint if met
            for constraint in strategy_json["structural_constraints"]:
                if constraint["type"] == "sequence" and constraint["details"]["target"] == "Silyl protective group" and constraint["details"]["start_position"] == "starting_material" and constraint["details"]["end_position"] == "late_stage":
                    if constraint not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint)
            # Add count constraint if met
            for constraint in strategy_json["structural_constraints"]:
                if constraint["type"] == "count" and constraint["details"]["target"] == "synthesis_stages_with_silyl_group" and constraint["details"]["operator"] == ">=" and constraint["details"]["value"] == 2:
                    if constraint not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint)
            print(f"Detected silyl protection maintained from starting material to late stage")

    # Case 3: Protection is observed at multiple stages including late stage
    elif len(protected_mol_at_depth) >= 2 and protection_maintained:
        depths = sorted(protected_mol_at_depth.keys())
        has_protection_strategy = True
        # Add count constraint if met
        for constraint in strategy_json["structural_constraints"]:
            if constraint["type"] == "count" and constraint["details"]["target"] == "synthesis_stages_with_silyl_group" and constraint["details"]["operator"] == ">=" and constraint["details"]["value"] == 2:
                if constraint not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(constraint)
        print(f"Detected silyl protection at multiple stages including late stage")

    return has_protection_strategy, findings_json
