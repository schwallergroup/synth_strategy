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


THIENOPYRIDINE_ISOMERS_SMARTS = [
    "c1csc2c1ncccc2",  # Thieno[2,3-b]pyridine
    "c1ccnc2c1scc2",  # Thieno[3,2-b]pyridine
    "c1cncc2c1scc2",  # Thieno[2,3-c]pyridine
    "c1cncc2c1csc2",  # Thieno[3,2-c]pyridine
]

def main(route) -> Tuple[bool, Dict]:
    """Detects if a synthesis maintains a specific thienopyridine core throughout its steps. The strategy is confirmed if the core is present in the final product and at least one preceding intermediate. The core's presence is verified by matching against a predefined list of thienopyridine isomer SMARTS patterns."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    core_present_in_final = False
    core_present_in_intermediates = False

    def has_thiophene_fused_pyridine(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        for smarts in THIENOPYRIDINE_ISOMERS_SMARTS:
            pattern = Chem.MolFromSmarts(smarts)
            if mol.HasSubstructMatch(pattern):
                return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal core_present_in_final, core_present_in_intermediates, findings_json

        if node.get("type") == "mol" and "smiles" in node:
            if has_thiophene_fused_pyridine(node["smiles"]):
                # Add to atomic checks for ring systems
                if "thienopyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("thienopyridine")

                if depth == 0:
                    core_present_in_final = True
                    # Add to structural constraints for last_stage
                    constraint_obj = {
                        "type": "positional",
                        "details": {
                            "target": "thienopyridine",
                            "position": "last_stage"
                        }
                    }
                    if constraint_obj not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint_obj)
                else:
                    core_present_in_intermediates = True
                    # Add to structural constraints for not_last_stage
                    constraint_obj = {
                        "type": "positional",
                        "details": {
                            "target": "thienopyridine",
                            "position": "not_last_stage"
                        }
                    }
                    if constraint_obj not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint_obj)

        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            if node.get("type") == "reaction":
                # If current node is a reaction, depth remains the same for children
                new_depth = depth
            else:
                # If current node is not a reaction (e.g., chemical), depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = core_present_in_final and core_present_in_intermediates
    return result, findings_json
