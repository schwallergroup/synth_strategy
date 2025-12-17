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


ISOCYANATE_TO_UREA_REACTIONS = [
    "Urea synthesis via isocyanate and primary amine",
    "Urea synthesis via isocyanate and secondary amine",
    "Urea synthesis via isocyanate and diazo",
    "Urea synthesis via isocyanate and sulfonamide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves a late-stage urea formation from an isocyanate.
    This is defined as a reaction matching a specific list of named urea syntheses occurring at depth 0 or 1.
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

    isocyanate_urea_at_late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal isocyanate_urea_at_late_stage, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check if this is an isocyanate-based urea formation reaction
            current_node_is_isocyanate_reaction = False
            for r in ISOCYANATE_TO_UREA_REACTIONS:
                if checker.check_reaction(r, rsmi):
                    current_node_is_isocyanate_reaction = True
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)

            # If this is at depth 0 or 1 (late stage) and involves isocyanate urea formation
            if depth <= 1 and current_node_is_isocyanate_reaction:
                isocyanate_urea_at_late_stage = True
                # Add the structural constraint if detected
                if {"type": "positional", "details": {"target": "Urea synthesis via isocyanate", "position": "late_stage (depth <= 1)"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Urea synthesis via isocyanate",
                            "position": "late_stage (depth <= 1)"
                        }
                    })

        # Traverse children (depth increases as we go backward in synthesis)
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            # Depth remains the same when traversing from reaction to chemical
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return isocyanate_urea_at_late_stage, findings_json
