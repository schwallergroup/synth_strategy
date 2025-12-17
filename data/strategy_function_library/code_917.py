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


SNAR_REACTION_TYPES = [
    "heteroaromatic_nuc_sub",
    "nucl_sub_aromatic_ortho_nitro",
    "nucl_sub_aromatic_para_nitro",
]

AMINE_NUCLEOPHILES = ["Aniline", "Primary amine", "Secondary amine"]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage SNAr reaction with an amine nucleophile.
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

    found = False
    is_snar_type_detected = False
    amine_nucleophile_in_reactant_detected = False
    late_stage_condition_met = False

    def dfs_traverse(node, depth=0):
        nonlocal found, is_snar_type_detected, amine_nucleophile_in_reactant_detected, late_stage_condition_met, findings_json

        if found:
            return

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if this is a late-stage reaction (depth <= 2)
            if depth <= 2:
                late_stage_condition_met = True
                # Check if this is a recognized SNAr reaction type
                current_is_snar_type = False
                for name in SNAR_REACTION_TYPES:
                    if checker.check_reaction(name, rsmi):
                        current_is_snar_type = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                
                if current_is_snar_type:
                    is_snar_type_detected = True
                    # Confirm an amine nucleophile was used
                    current_amine_found = False
                    for fg in AMINE_NUCLEOPHILES:
                        for r in reactants:
                            if checker.check_fg(fg, r):
                                current_amine_found = True
                                if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(fg)
                                break
                        if current_amine_found:
                            break

                    if current_amine_found:
                        amine_nucleophile_in_reactant_detected = True
                        # If both conditions are met, record the co-occurrence constraint
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "any_SNAR_reaction",
                                    "amine_nucleophile_in_reactant"
                                ],
                                "scope": "same_reaction_step",
                                "description": "A reaction must be one of the specified SNAr types and one of its reactants must contain one of the specified amine functional groups."
                            }
                        })
                        found = True
                        return

        # Continue DFS traversal
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if found and late_stage_condition_met:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "co-occurring_SNAr_and_amine_reaction",
                "position": "late_stage",
                "condition": "depth <= 2"
            }
        })

    return found, findings_json
