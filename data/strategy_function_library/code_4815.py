from typing import Tuple, Dict, List
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


N_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Phthalimide deprotection",
    "N-glutarimide deprotection",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route contains a late-stage N-deprotection step.
    Late stage is defined as occurring in the first half of the synthesis (low depth).
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

    has_late_deprotection = False
    max_depth = 0

    # First, find the maximum depth to determine what "late stage" means
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    find_max_depth(route)
    print(f"Maximum depth of synthesis: {max_depth}")

    # Now check for deprotection in the late stage
    def dfs_traverse(node, depth=0):
        nonlocal has_late_deprotection, findings_json

        # Consider late stage as the first half of the synthesis (lower depth numbers)
        is_late_stage = depth <= max_depth / 2

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check for specific N-deprotection reactions
            is_deprotection = False
            for reaction_type in N_DEPROTECTION_REACTIONS:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Found {reaction_type} reaction at depth {depth}")
                    is_deprotection = True
                    if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)

                    # Check if this is a late-stage deprotection
                    if is_late_stage:
                        print(f"This is a late-stage deprotection (depth {depth} <= {max_depth/2})")
                        has_late_deprotection = True
                        # Add structural constraint for positional detection
                        if {"type": "positional", "details": {"target": N_DEPROTECTION_REACTIONS, "position": "first_half"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": N_DEPROTECTION_REACTIONS, "position": "first_half"}})
                    else:
                        print(
                            f"This is an early-stage deprotection (depth {depth} > {max_depth/2})"
                        )

                    break

            # If no specific reaction type matched, check for protected groups in reactants
            # and their absence in products
            if not is_deprotection and is_late_stage:
                # Check for Boc or Phthalimide groups in reactants
                for reactant in reactants:
                    # Check for Boc group (might be missed by the checker)
                    has_boc = checker.check_fg("Boc", reactant)
                    has_phthalimide = checker.check_fg("Phthalimide", reactant)

                    if has_boc:
                        if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Boc")
                    if has_phthalimide:
                        if "Phthalimide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Phthalimide")

                    if has_boc or has_phthalimide:
                        protecting_group = "Boc" if has_boc else "Phthalimide"
                        print(f"Found protected group {protecting_group} in reactant: {reactant}")

                        # Check if the product no longer has this protection
                        product_has_boc = checker.check_fg(
                            "Boc", product
                        )
                        product_has_phthalimide = checker.check_fg("Phthalimide", product)

                        if (has_boc and not product_has_boc) or (
                            has_phthalimide and not product_has_phthalimide
                        ):
                            print(f"Product no longer has {protecting_group} protection: {product}")

                            # Verify it's an amine deprotection by checking for amine in product
                            amine_types = [
                                "Primary amine",
                                "Secondary amine",
                                "Tertiary amine",
                                "Aniline",
                            ]
                            found_amine_type = None
                            for amine_type in amine_types:
                                if checker.check_fg(amine_type, product):
                                    if amine_type not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append(amine_type)
                                    found_amine_type = amine_type
                                    break

                            if found_amine_type:
                                print(f"Confirmed N-deprotection at depth {depth} (late stage)")
                                has_late_deprotection = True
                                # Add structural constraint for co-occurrence
                                if protecting_group == "Boc":
                                    constraint = {"type": "co-occurrence", "details": {"targets": ["Boc", "Primary amine", "Secondary amine", "Tertiary amine", "Aniline"], "position": "first_half"}}
                                elif protecting_group == "Phthalimide":
                                    constraint = {"type": "co-occurrence", "details": {"targets": ["Phthalimide", "Primary amine", "Secondary amine", "Tertiary amine", "Aniline"], "position": "first_half"}}
                                if constraint not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append(constraint)
                                break

                # Also check for nitro reduction to amine (another form of "deprotection")
                if not has_late_deprotection:
                    for reactant in reactants:
                        if checker.check_fg("Nitro group", reactant):
                            if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

                        if checker.check_fg("Nitro group", reactant) and not checker.check_fg(
                            "Nitro group", product
                        ):
                            amine_types = [
                                "Primary amine",
                                "Secondary amine",
                                "Tertiary amine",
                                "Aniline",
                            ]
                            found_amine_type = None
                            for amine_type in amine_types:
                                if checker.check_fg(amine_type, product):
                                    if amine_type not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append(amine_type)
                                    found_amine_type = amine_type
                                    break

                            if found_amine_type:
                                print(
                                    f"Found nitro reduction to amine at depth {depth} (late stage)"
                                )
                                has_late_deprotection = True
                                # Add structural constraint for co-occurrence
                                constraint = {"type": "co-occurrence", "details": {"targets": ["Nitro group", "Primary amine", "Secondary amine", "Tertiary amine", "Aniline"], "position": "first_half"}}
                                if constraint not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append(constraint)
                                break

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_late_deprotection, findings_json
