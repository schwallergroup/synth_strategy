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


SULFONAMIDE_FORMATION_REACTIONS = [
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage sulfonamide formation, defined as occurring within the final 4 steps. A reaction is flagged if a sulfonamide group is formed from reactants lacking one, and the transformation is identified as a specific named reaction (from `SULFONAMIDE_FORMATION_REACTIONS`) or by the presence of sulfonyl halide and amine reactants.
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

    sulfonamide_formation_detected = False
    depth_of_formation = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_formation_detected, depth_of_formation, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                if checker.check_fg("Sulfonamide", product):
                    findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")

                    reactants_with_sulfonamide = [
                        r for r in reactants if checker.check_fg("Sulfonamide", r)
                    ]

                    if not reactants_with_sulfonamide:
                        # This corresponds to the 'negation' structural constraint
                        findings_json["structural_constraints"].append({
                            "type": "negation",
                            "details": {
                                "target": "Sulfonamide in reactants"
                            }
                        })

                        is_sulfonamide_reaction = False
                        for name in SULFONAMIDE_FORMATION_REACTIONS:
                            if checker.check_reaction(name, rsmi):
                                is_sulfonamide_reaction = True
                                findings_json["atomic_checks"]["named_reactions"].append(name)
                                break

                        has_sulfonyl_chloride = False
                        for r in reactants:
                            if checker.check_fg("Sulfonyl halide", r):
                                has_sulfonyl_chloride = True
                                findings_json["atomic_checks"]["functional_groups"].append("Sulfonyl halide")
                                break

                        has_amine = False
                        for r in reactants:
                            if checker.check_fg("Primary amine", r):
                                has_amine = True
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                            if checker.check_fg("Secondary amine", r):
                                has_amine = True
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                            if has_amine: # Optimization: if either amine is found, no need to check further
                                break

                        if is_sulfonamide_reaction or (has_sulfonyl_chloride and has_amine):
                            sulfonamide_formation_detected = True
                            depth_of_formation = min(depth_of_formation, depth)
            except Exception:
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = sulfonamide_formation_detected and depth_of_formation <= 4
    if result:
        # This corresponds to the 'positional' structural constraint
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "sulfonamide_formation",
                "position": "within_final_4_steps"
            }
        })

    return result, findings_json
