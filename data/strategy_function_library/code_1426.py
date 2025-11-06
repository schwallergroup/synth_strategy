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


# Refactored lists for enumeration
BOC_PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
]

BOC_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a Boc protecting group is maintained throughout most of the synthesis.
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

    boc_presence_by_depth = {}
    max_depth = 0
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_presence_by_depth, max_depth, findings_json
        max_depth = max(max_depth, depth)

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol_smiles = node["smiles"]
                # Check if Boc group is present
                has_boc = checker.check_fg("Boc", mol_smiles)
                if has_boc:
                    if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Boc")
                print(
                    f"Depth {depth}: {'Boc found in' if has_boc else 'No Boc found in'} molecule {mol_smiles}"
                )
                boc_presence_by_depth[depth] = has_boc

            except Exception as e:
                print(f"Error processing molecule at depth {depth}: {e}")

        # Check reaction nodes to verify Boc is maintained through reactions
        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if Boc is present in reactants and products
                reactants_have_boc = any(checker.check_fg("Boc", r) for r in reactants)
                product_has_boc = checker.check_fg("Boc", product)

                # Check if this is a Boc protection or deprotection reaction
                is_boc_protection = False
                for r_name in BOC_PROTECTION_REACTIONS:
                    if checker.check_reaction(r_name, rsmi):
                        is_boc_protection = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        break

                is_boc_deprotection = False
                for r_name in BOC_DEPROTEPTION_REACTIONS:
                    if checker.check_reaction(r_name, rsmi):
                        is_boc_deprotection = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        # Record deprotection as a negative finding
                        if {"type": "negation", "details": {"description": "The route should not contain a Boc deprotection reaction, which would break the maintained protection.", "targets": [r_name]}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "negation", "details": {"description": "The route should not contain a Boc deprotection reaction, which would break the maintained protection.", "targets": [r_name]}})
                        break

                if is_boc_protection:
                    print(f"Depth {depth}: Boc protection reaction detected")
                    boc_presence_by_depth[depth] = True
                elif is_boc_deprotection:
                    print(
                        f"Depth {depth}: Boc deprotection reaction detected - this breaks the maintained protection"
                    )
                    boc_presence_by_depth[depth] = False
                else:
                    # For other reactions, check if Boc is maintained
                    if reactants_have_boc and not product_has_boc:
                        print(f"Depth {depth}: Boc lost in reaction")
                        boc_presence_by_depth[depth] = False
                    elif not reactants_have_boc and product_has_boc:
                        print(f"Depth {depth}: Boc introduced in reaction")
                        boc_presence_by_depth[depth] = True
                    else:
                        print(
                            f"Depth {depth}: Boc {'maintained' if product_has_boc else 'not present'} in reaction"
                        )
                        boc_presence_by_depth[depth] = product_has_boc

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "mol"
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if Boc is present in at least 70% of the synthesis steps
    # Give more weight to late-stage steps (lower depth values)
    boc_steps = 0
    total_steps = 0

    for depth, present in boc_presence_by_depth.items():
        # Late-stage steps have lower depth values
        # Give full weight to final product and immediate precursors, half weight to earlier steps
        weight = 1.0 if depth <= 2 else 0.5

        if present:
            boc_steps += weight
        total_steps += weight

    # Calculate the weighted percentage
    if total_steps > 0:
        boc_percentage = boc_steps / total_steps
        print(f"Boc protection maintained with weighted percentage: {boc_percentage:.2f}")

        # Check if Boc is present in at least 70% of steps (weighted)
        if boc_percentage >= 0.7:
            # Additional check: verify Boc is present in the final product (depth 0)
            if 0 in boc_presence_by_depth and boc_presence_by_depth[0]:
                print("Boc protection maintained throughout synthesis and present in final product")
                result = True
                if {"type": "positional", "details": {"target": "Boc", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Boc", "position": "last_stage"}})
                if {"type": "count", "details": {"target": "weighted_boc_presence_ratio", "operator": ">=", "value": 0.7}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "count", "details": {"target": "weighted_boc_presence_ratio", "operator": ">=", "value": 0.7}})
            else:
                print("Boc not present in final product, so not maintained throughout")
                result = False
        else:
            print(f"Boc protection not maintained in enough steps ({boc_percentage:.2f} < 0.7)")
            result = False

    else:
        print("No steps with Boc protection found")
        result = False

    return result, findings_json
