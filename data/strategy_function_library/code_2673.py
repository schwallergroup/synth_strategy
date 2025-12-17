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


SULFONAMIDE_FORMATION_REACTIONS = [
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
    "sulfon_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving late-stage formation of a sulfonamide group.
    Late-stage means in the second half of the synthesis (higher depth values in retrosynthesis).
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

    found_sulfonamide_formation = False
    sulfonamide_formation_depth = float("inf")  # Initialize to infinity
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_sulfonamide_formation, sulfonamide_formation_depth, max_depth, findings_json

        # Track maximum depth to determine synthesis length
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")

            # Check if product has sulfonamide
            product_has_sulfonamide = checker.check_fg("Sulfonamide", product_part)
            if product_has_sulfonamide and "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")

            # Check if this is a sulfonamide formation reaction
            is_sulfonamide_reaction = False
            for r in SULFONAMIDE_FORMATION_REACTIONS:
                if checker.check_reaction(r, rsmi):
                    is_sulfonamide_reaction = True
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)
            
            if is_sulfonamide_reaction and product_has_sulfonamide:
                # Check if any reactant doesn't have the sulfonamide
                has_reactant_without_sulfonamide = any(
                    not checker.check_fg("Sulfonamide", reactant) for reactant in reactants
                )

                if has_reactant_without_sulfonamide:
                    # If we find a sulfonamide formation at a lower depth (later stage),
                    # or if this is the first one we've found, record it
                    if not found_sulfonamide_formation or depth < sulfonamide_formation_depth:
                        found_sulfonamide_formation = True
                        sulfonamide_formation_depth = depth
                        # print(f"Found sulfonamide formation at depth {depth}, reaction: {rsmi}")

            # Additional check for sulfonamide formation by looking at functional groups
            elif product_has_sulfonamide:
                # Check if any reactant doesn't have the sulfonamide
                has_reactant_without_sulfonamide = any(
                    not checker.check_fg("Sulfonamide", reactant) for reactant in reactants
                )

                if has_reactant_without_sulfonamide:
                    # Check if reactants have sulfonyl chloride and amine
                    has_sulfonyl_chloride = False
                    for reactant in reactants:
                        if checker.check_fg("Sulfonyl halide", reactant):
                            has_sulfonyl_chloride = True
                            if "Sulfonyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Sulfonyl halide")
                            break

                    has_amine = False
                    for reactant in reactants:
                        if checker.check_fg("Primary amine", reactant):
                            has_amine = True
                            if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        elif checker.check_fg("Secondary amine", reactant):
                            has_amine = True
                            if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                        elif checker.check_fg("Aniline", reactant):
                            has_amine = True
                            if "Aniline" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Aniline")
                        if has_amine: break

                    if has_sulfonyl_chloride and has_amine:
                        # This is likely a sulfonamide formation reaction not captured by the reaction checks
                        if not found_sulfonamide_formation or depth < sulfonamide_formation_depth:
                            found_sulfonamide_formation = True
                            sulfonamide_formation_depth = depth
                            if "sulfonamide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("sulfonamide_formation")
                            # print(
                            #     f"Found sulfonamide formation via functional groups at depth {depth}"
                            # )

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                # From reaction to chemical, depth remains the same
                dfs_traverse(child, depth)
            else:
                # From chemical to reaction, depth increases
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if sulfonamide formation is late-stage (in the second half of synthesis)
    # In retrosynthesis, higher depth values are later in the forward synthesis
    late_stage = False
    if found_sulfonamide_formation and max_depth > 0:
        midpoint = max_depth / 2
        # Late stage means depth is greater than or equal to midpoint (deeper in the tree)
        late_stage = sulfonamide_formation_depth >= midpoint
        # print(
        #     f"Max depth: {max_depth}, Midpoint: {midpoint}, Sulfonamide formation depth: {sulfonamide_formation_depth}"
        # )
        # print(
        #     f"Is late stage? {late_stage} (depth {sulfonamide_formation_depth} >= midpoint {midpoint})"
        # )
        if late_stage:
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "sulfonamide_formation",
                    "position": "second_half"
                }
            })
    # else:
        # print(f"Max depth: {max_depth}, Sulfonamide formation not found or max_depth is 0")

    return late_stage, findings_json