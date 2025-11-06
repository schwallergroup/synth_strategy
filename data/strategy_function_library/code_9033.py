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
    This function detects if the synthetic route employs a late-stage nitro reduction strategy,
    where a nitro group is carried through most of the synthesis and reduced to an amine
    only in the final step (depth 1).
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

    # Flag to track if we found a nitro reduction at depth 1
    found_late_stage_nitro_reduction = False
    # Track if we found nitro groups in earlier steps
    found_nitro_in_earlier_steps = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_nitro_reduction, found_nitro_in_earlier_steps, findings_json

        # Get depth from metadata if available
        if "depth" in node.get("metadata", {}):
            depth = node["metadata"]["depth"]

        print(f"Traversing node type: {node.get('type')} at depth: {depth}")

        # Check if this is a molecule node
        if node.get("type") == "mol":
            mol_smiles = node.get("smiles", "")

            # Check for nitro groups in earlier steps (depth > 1)
            if depth > 1 and checker.check_fg("Nitro group", mol_smiles):
                print(f"Found nitro group at depth {depth}: {mol_smiles}")
                found_nitro_in_earlier_steps = True
                if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
                # Add structural constraint for nitro group not in last stage
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Nitro group", "position": "not_last_stage"}})

        # Check if this is a reaction node
        elif node.get("type") == "reaction":
            # Check if this is the final reaction (depth 1)
            if depth == 1:
                # Extract reactants and product
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if rsmi:
                    try:
                        parts = rsmi.split(">")
                        if len(parts) >= 3:
                            reactants_smiles = parts[0].split(".")
                            product_smiles = parts[2]

                            print(f"Analyzing reaction at depth 1: {rsmi}")

                            # Check for nitro groups in reactants
                            nitro_reactants = [
                                r for r in reactants_smiles if checker.check_fg("Nitro group", r)
                            ]
                            if nitro_reactants:
                                if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

                            # Check for amine groups in product
                            has_amine_product = False
                            amine_types = ["Primary amine", "Secondary amine", "Tertiary amine", "Aniline"]
                            for amine_type in amine_types:
                                if checker.check_fg(amine_type, product_smiles):
                                    has_amine_product = True
                                    if amine_type not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append(amine_type)

                            # Check if this is a nitro reduction reaction
                            is_nitro_reduction = checker.check_reaction(
                                "Reduction of nitro groups to amines", rsmi
                            )
                            if is_nitro_reduction:
                                if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

                            # If not detected by reaction type, check by functional group changes
                            if not is_nitro_reduction and nitro_reactants and has_amine_product:
                                # Check if nitro group is gone in product
                                if not checker.check_fg("Nitro group", product_smiles):
                                    print(
                                        "Detected nitro reduction based on functional group changes"
                                    )
                                    is_nitro_reduction = True
                                    # Even if not by name, it's a reduction of nitro groups to amines by FG change
                                    if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                                        findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

                            if is_nitro_reduction and nitro_reactants:
                                print("Confirmed nitro reduction at final step")
                                found_late_stage_nitro_reduction = True
                                # Add structural constraint for reduction at last stage
                                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Reduction of nitro groups to amines", "position": "last_stage"}})
                                # Add structural constraint for negation of nitro group in product
                                if not checker.check_fg("Nitro group", product_smiles):
                                    findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "Nitro group", "scope": "product_of_last_stage_reaction"}})

                    except Exception as e:
                        print(f"Error processing reaction SMILES: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node.get("type") != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Continue traversing children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True only if we found both a late-stage nitro reduction and nitro groups in earlier steps
    result = found_late_stage_nitro_reduction and found_nitro_in_earlier_steps

    if found_late_stage_nitro_reduction and found_nitro_in_earlier_steps:
        # Add co-occurrence constraint if both conditions are met
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Reduction of nitro groups to amines", "Nitro group"]}})

    print(
        f"Late stage nitro reduction: {found_late_stage_nitro_reduction}, Nitro in earlier steps: {found_nitro_in_earlier_steps}"
    )
    return result, findings_json
