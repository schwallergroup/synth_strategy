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
    Detects late-stage amide formation (depth=1) via acylation of a nitrogen nucleophile by an acyl halide or related species.
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

    late_amide_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal late_amide_formation_detected, findings_json

        if node["type"] == "reaction" and depth <= 1:
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["mapped_reaction_smiles"]

                    reaction_name = "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N"
                    is_acylation = checker.check_reaction(
                        reaction_name,
                        rsmi,
                    )

                    if is_acylation:
                        late_amide_formation_detected = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        # Add structural constraint if detected at late stage (depth <= 1)
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": reaction_name,
                                "position": "late_stage"
                            }
                        })
                        return
            except Exception:
                pass

        for child in node.get("children", []):
            if not late_amide_formation_detected:
                # New logic for depth calculation
                new_depth = depth
                if node["type"] != "reaction": # If current node is chemical, depth increases
                    new_depth = depth + 1
                dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return late_amide_formation_detected, findings_json