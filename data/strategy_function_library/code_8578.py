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


SULFONAMIDE_NAMED_REACTIONS = [
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage sulfonamide formation, defined as occurring within the final three synthetic steps. 
    It first attempts to identify the transformation by checking for specific named reactions from the `SULFONAMIDE_NAMED_REACTIONS` list. 
    If no named reaction is matched, it falls back to a general analysis, confirming the consumption of an amine and a sulfonyl halide to form a new sulfonamide group that was not present in the reactants.
    """
    print(f"Analyzing route for late-stage sulfonamide formation")
    late_stage_sulfonamide = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_sulfonamide, findings_json

        if node["type"] == "reaction" and depth <= 3:
            # Record positional constraint if this reaction is within the last three stages
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "sulfonamide_formation",
                    "position": "last_three_stages"
                }
            })

            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is a sulfonamide formation reaction - try specific reactions first
                is_sulfonamide_formation = False
                for rxn in SULFONAMIDE_NAMED_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        is_sulfonamide_formation = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break

                print(
                    f"Named reaction check at depth {depth}: found={is_sulfonamide_formation}"
                )

                # If specific reaction checks fail, look for general pattern
                if not is_sulfonamide_formation:
                    has_sulfonyl_chloride = False
                    for reactant in reactants:
                        if checker.check_fg("Sulfonyl halide", reactant):
                            has_sulfonyl_chloride = True
                            findings_json["atomic_checks"]["functional_groups"].append("Sulfonyl halide")
                            break

                    has_amine = False
                    for reactant in reactants:
                        if checker.check_fg("Primary amine", reactant):
                            has_amine = True
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        if checker.check_fg("Secondary amine", reactant):
                            has_amine = True
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                        if has_amine: # If either primary or secondary amine is found, we can break
                            break

                    has_sulfonamide_product = checker.check_fg("Sulfonamide", product)
                    if has_sulfonamide_product:
                        findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")

                    sulfonamide_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg("Sulfonamide", reactant):
                            sulfonamide_in_reactants = True
                            # Do not add to findings_json as this is a negation check
                            break

                    if (
                        has_sulfonyl_chloride
                        and has_amine
                        and has_sulfonamide_product
                        and not sulfonamide_in_reactants
                    ):
                        print(f"Detected sulfonamide formation by component analysis at depth {depth}")
                        is_sulfonamide_formation = True
                        # Record co-occurrence constraint
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "Sulfonyl halide_in_reactant",
                                    "Amine_in_reactant",
                                    "Sulfonamide_in_product"
                                ],
                                "scope": "reaction"
                            }
                        })
                        # Record negation constraint
                        findings_json["structural_constraints"].append({
                            "type": "negation",
                            "details": {
                                "target": "Sulfonamide_in_reactant",
                                "scope": "reaction"
                            }
                        })

                # Confirm if a sulfonamide formation was found by either method
                if is_sulfonamide_formation:
                    print(f"Confirmed late-stage sulfonamide formation at depth {depth}")
                    late_stage_sulfonamide = True

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)
    return late_stage_sulfonamide, findings_json
