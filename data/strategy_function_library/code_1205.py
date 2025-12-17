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


AMINE_TYPES_FOR_SULFONAMIDE = ["Primary amine", "Secondary amine", "Aniline"]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route involves late-stage sulfonamide formation.
    Looks for sulfonamide formation in the second half of the synthesis.
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

    result = False
    max_depth = 0
    sulfonamide_depth = None

    # First pass to find max depth
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    # Second pass to find sulfonamide formation
    def find_sulfonamide(node, depth=0):
        nonlocal sulfonamide_depth, findings_json, result

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]
            reactants = reactants_str.split(".")

            # Check for key functional groups
            has_sulfonyl_halide = any(checker.check_fg("Sulfonyl halide", r) for r in reactants)
            if has_sulfonyl_halide:
                if "Sulfonyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Sulfonyl halide")

            has_sulfonamide = checker.check_fg("Sulfonamide", product_str)
            if has_sulfonamide:
                if "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")

            no_sulfonamide_in_reactants = not any(
                checker.check_fg("Sulfonamide", r) for r in reactants
            )
            if not no_sulfonamide_in_reactants:
                # This means sulfonamide was found in reactants, which is a negation constraint violation
                # We record the negation constraint if it's violated, as it's part of the strategy.
                if {"type": "negation", "details": {"target": "Sulfonamide", "context": "in reactants of formation step"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "Sulfonamide", "context": "in reactants of formation step"}})

            has_amine = False
            for amine_type in AMINE_TYPES_FOR_SULFONAMIDE:
                if any(checker.check_fg(amine_type, r) for r in reactants):
                    has_amine = True
                    if amine_type not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(amine_type)

            # Determine if this is a sulfonamide formation
            is_forming_sulfonamide = (
                has_sulfonamide
                and has_sulfonyl_halide
                and has_amine
                and no_sulfonamide_in_reactants
            )

            if is_forming_sulfonamide:
                if "Sulfonamide formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Sulfonamide formation")
                # Track the earliest (lowest depth) sulfonamide formation
                if sulfonamide_depth is None or depth < sulfonamide_depth:
                    sulfonamide_depth = depth

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is reaction, depth remains same for children (chemical nodes)
                find_sulfonamide(child, depth)
            else:
                # If current node is chemical, depth increases for children (reaction nodes)
                find_sulfonamide(child, depth + 1)

    find_max_depth(route)
    find_sulfonamide(route)

    # Check if sulfonamide formation occurs in second half of synthesis
    if sulfonamide_depth is not None:
        is_late_stage = sulfonamide_depth <= (max_depth // 2)
        if is_late_stage:
            result = True
            if {"type": "positional", "details": {"target": "Sulfonamide formation", "position": "second_half"}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Sulfonamide formation", "position": "second_half"}})

    return result, findings_json
