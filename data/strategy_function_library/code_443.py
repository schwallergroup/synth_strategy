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
    This function detects a synthetic strategy involving amine to nitrile conversion
    in the early stages of synthesis.
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

    amine_to_nitrile_found = False
    # Define what constitutes "early stages" of synthesis
    EARLY_STAGE_THRESHOLD = 3

    def dfs_traverse(node, depth=0):
        nonlocal amine_to_nitrile_found, findings_json

        if node["type"] == "reaction" and depth < EARLY_STAGE_THRESHOLD:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                amine_in_reactants = False
                for reactant in reactants:
                    if checker.check_fg("Primary amine", reactant):
                        amine_in_reactants = True
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")

                nitrile_in_product = checker.check_fg("Nitrile", product)
                if nitrile_in_product:
                    if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

                # Forward direction: amine to nitrile
                if amine_in_reactants and nitrile_in_product:
                    reaction_found_in_pathway = False
                    for rxn_type in ["Amine to azide", "Azide to nitrile"]:
                        if checker.check_reaction(rxn_type, rsmi):
                            reaction_found_in_pathway = True
                            if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_type)

                    if reaction_found_in_pathway:
                        amine_to_nitrile_found = True
                        # Record the co-occurrence structural constraint if all conditions met
                        if amine_in_reactants and nitrile_in_product and reaction_found_in_pathway:
                            co_occurrence_constraint = {
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "Primary amine",
                                        "Nitrile",
                                        "Amine to azide",
                                        "Azide to nitrile"
                                    ],
                                    "scope": "reaction_step",
                                    "logic": "(target_1 AND target_2) AND (target_3 OR target_4)"
                                }
                            }
                            if co_occurrence_constraint not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(co_occurrence_constraint)

        # Continue traversing children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Record positional constraint if amine_to_nitrile_found is True and depth is within early stage
    # This positional constraint is based on the overall finding, not a specific node.
    if amine_to_nitrile_found:
        positional_constraint = {
            "type": "positional",
            "details": {
                "target": "amine_to_nitrile_conversion",
                "position": "late_stage",
                "max_depth": 3
            }
        }
        # The original logic checks for 'early stages' (depth < EARLY_STAGE_THRESHOLD)
        # The JSON constraint says 'late_stage' with max_depth 3. This implies if found within depth 3, it's 'late_stage' for this constraint.
        # Assuming the 'late_stage' here means 'found within the first 3 steps' as per the original function's intent.
        # If the strategy is found, and it's within the early stage threshold, then this constraint is met.
        # The original function's `EARLY_STAGE_THRESHOLD` is 3, and the JSON `max_depth` is 3.
        # So, if `amine_to_nitrile_found` is true, it means it was found at depth < 3. This matches the positional constraint.
        if positional_constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(positional_constraint)

    return amine_to_nitrile_found, findings_json
