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
    This function detects if the route includes amide formation following a nitro reduction.

    The process involves:
    1. Reduction of a nitro group to an amine
    2. A subsequent or concurrent amide formation step
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

    # Track both nitro reduction and subsequent lactam formation
    nitro_reduction_nodes = []
    lactam_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal lactam_formation_found, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            try:
                # Extract reactants and product using the correct format
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Check for nitro reduction
                has_nitro_reactant = checker.check_fg("Nitro group", reactants_part)
                if has_nitro_reactant:
                    if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

                has_amine_product = (
                    checker.check_fg("Primary amine", product_part)
                    or checker.check_fg("Secondary amine", product_part)
                )
                if checker.check_fg("Primary amine", product_part) and "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                if checker.check_fg("Secondary amine", product_part) and "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

                nitro_reduced = has_nitro_reactant and not checker.check_fg(
                    "Nitro group", product_part
                )
                is_reduction_reaction = checker.check_reaction(
                    "Reduction of nitro groups to amines", rsmi
                )
                if is_reduction_reaction:
                    if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

                if (
                    has_nitro_reactant
                    and has_amine_product
                    and nitro_reduced
                    and is_reduction_reaction
                ):
                    print(f"Detected nitro reduction at depth {depth}: {rsmi}")
                    nitro_reduction_nodes.append((node, depth))

                # Check for amide formation
                has_amide_product = (
                    checker.check_fg("Primary amide", product_part)
                    or checker.check_fg("Secondary amide", product_part)
                    or checker.check_fg("Tertiary amide", product_part)
                )
                if checker.check_fg("Primary amide", product_part) and "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                if checker.check_fg("Secondary amide", product_part) and "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                if checker.check_fg("Tertiary amide", product_part) and "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                # Check for direct amide formation via nitro reduction in a single step
                if (
                    has_nitro_reactant
                    and has_amide_product
                    and nitro_reduced
                ):
                    print(
                        f"Detected direct lactam formation via nitro reduction at depth {depth}: {rsmi}"
                    )
                    lactam_formation_found = True
                    if {"type": "sequence", "details": {"ordered_events": ["Reduction of nitro groups to amines", "amide_formation"], "description": "An amide formation event, identified by the appearance of an amide functional group, must occur in the same reaction step as a nitro reduction or in a subsequent step."}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["Reduction of nitro groups to amines", "amide_formation"], "description": "An amide formation event, identified by the appearance of an amide functional group, must occur in the same reaction step as a nitro reduction or in a subsequent step."}})

                # Check for amide formation in a step after nitro reduction
                elif has_amide_product and not has_nitro_reactant:
                    # Check if this step follows a nitro reduction step
                    for nitro_node, nitro_depth in nitro_reduction_nodes:
                        # If this step is later in the synthesis (lower depth) than a nitro reduction
                        if depth < nitro_depth:
                            print(
                                f"Detected lactam formation following nitro reduction at depth {depth}: {rsmi}"
                            )
                            print(f"  Related to nitro reduction at depth {nitro_depth}")
                            lactam_formation_found = True
                            if {"type": "sequence", "details": {"ordered_events": ["Reduction of nitro groups to amines", "amide_formation"], "description": "An amide formation event, identified by the appearance of an amide functional group, must occur in the same reaction step as a nitro reduction or in a subsequent step."}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["Reduction of nitro groups to amines", "amide_formation"], "description": "An amide formation event, identified by the appearance of an amide functional group, must occur in the same reaction step as a nitro reduction or in a subsequent step."}})
                            break

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return lactam_formation_found, findings_json
