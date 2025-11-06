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


AROMATIC_HALOGENATION_REACTIONS = [
    "Aromatic fluorination",
    "Aromatic chlorination",
    "Aromatic bromination",
    "Aromatic iodination",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves a late-stage aromatic halogenation, specifically checking for the following named reactions: Aromatic fluorination, Aromatic chlorination, Aromatic bromination, and Aromatic iodination.
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

    halogenation_detected = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal halogenation_detected, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is an aromatic halogenation reaction
                is_aromatic_halogenation = False
                for rxn_name in AROMATIC_HALOGENATION_REACTIONS:
                    if checker.check_reaction(rxn_name, rsmi):
                        is_aromatic_halogenation = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_name)

                # Determine if this is a late-stage reaction (in the final 1/3 of the synthesis)
                is_late_stage = depth <= max(1, max_depth // 3)

                if is_aromatic_halogenation and is_late_stage:
                    halogenation_detected = True
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": [
                                "Aromatic fluorination",
                                "Aromatic chlorination",
                                "Aromatic bromination",
                                "Aromatic iodination"
                            ],
                            "position": "late_stage"
                        }
                    })
                    print(
                        f"Late-stage aromatic halogenation detected at depth {depth} of {max_depth}"
                    )

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = False
    # Only consider it late-stage if the synthesis has at least 2 steps
    if max_depth >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "total_synthesis_steps",
                "operator": ">=",
                "value": 2
            }
        })
        if halogenation_detected:
            result = True

    return result, findings_json
