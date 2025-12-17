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


ESTERIFICATION_REACTIONS_OF_INTEREST = [
    "Esterification of Carboxylic Acids",
    "Schotten-Baumann to ester",
    "Transesterification",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route involves a late-stage esterification, based on a specific list of reaction types.
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

    late_stage_esterification = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_esterification, findings_json

        if (
            not late_stage_esterification and node["type"] == "reaction" and depth <= 1
        ):  # Only consider late-stage reactions (low depth)
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                for reaction_name in ESTERIFICATION_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(reaction_name, rsmi):
                        late_stage_esterification = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": [
                                    "Esterification of Carboxylic Acids",
                                    "Schotten-Baumann to ester",
                                    "Transesterification"
                                ],
                                "position": "depth <= 1"
                            }
                        })
                        break
            except Exception:
                pass

        # Continue traversing with modified depth
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return late_stage_esterification, findings_json
