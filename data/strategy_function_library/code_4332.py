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

def main(route) -> Tuple[bool, Dict]:
    """Checks for early or mid-synthesis Boc deprotection. The strategy is considered valid if the route contains no amines, or if any 'Boc amine deprotection' reactions occur only in the final step of the synthesis."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    has_amine_groups = False
    late_stage_deprotection_only = True

    max_depth = 0

    def get_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)
        for child in node.get("children", []):
            get_depth(child, current_depth + 1)

    get_depth(route)

    def dfs_traverse(node, depth=0, branch_path=None):
        nonlocal late_stage_deprotection_only, has_amine_groups, findings_json

        if branch_path is None:
            branch_path = []

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            has_primary_amine = checker.check_fg("Primary amine", mol_smiles)
            has_secondary_amine = checker.check_fg("Secondary amine", mol_smiles)

            if has_primary_amine:
                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
            if has_secondary_amine:
                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

            if has_primary_amine or has_secondary_amine:
                has_amine_groups = True

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rxn_smiles = node["metadata"]["mapped_reaction_smiles"]
            is_boc_deprotection = checker.check_reaction("Boc amine deprotection", rxn_smiles)

            if is_boc_deprotection:
                findings_json["atomic_checks"]["named_reactions"].append("Boc amine deprotection")
                # A deprotection is only allowed in the final step (depth=1).
                # Any deprotection at a greater depth is too early.
                if depth > 1:
                    late_stage_deprotection_only = False

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                dfs_traverse(child, depth, branch_path + [node])
            else: # Assuming 'mol' or other types that should increase depth
                dfs_traverse(child, depth + 1, branch_path + [node])

    dfs_traverse(route)

    result = False
    if not has_amine_groups:
        result = True
    else:
        result = late_stage_deprotection_only

    if result and has_amine_groups and late_stage_deprotection_only:
        # This condition implies that Boc deprotection occurred, but only in the last stage.
        # This corresponds to the structural constraint.
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Boc amine deprotection",
                "position": "last_stage"
            }
        })

    return result, findings_json
