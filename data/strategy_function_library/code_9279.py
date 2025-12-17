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
r_name_azide_nitrile = "Azide-nitrile click cycloaddition to tetrazole"

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
    Detects if the synthesis route involves late-stage tetrazole formation via an azide-nitrile cycloaddition. Late-stage is defined as occurring in the final two steps (depth <= 1).
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

    tetrazole_formation_detected = False
    tetrazole_formation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal tetrazole_formation_detected, tetrazole_formation_depth, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Direct check for the specific reaction type
                if checker.check_reaction(r_name_azide_nitrile, rsmi):
                    print(f"Detected tetrazole formation reaction at depth {depth}")
                    tetrazole_formation_detected = True
                    tetrazole_formation_depth = depth
                    if r_name_azide_nitrile not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name_azide_nitrile)
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if tetrazole formation is late-stage (depth 0 or 1)
    is_late_stage = (
        tetrazole_formation_detected
        and tetrazole_formation_depth is not None
        and tetrazole_formation_depth <= 1
    )

    if is_late_stage:
        print(f"Late-stage tetrazole formation detected at depth {tetrazole_formation_depth}")
        # Add the structural constraint if late-stage condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Azide-nitrile click cycloaddition to tetrazole",
                "position_type": "max_depth_from_product",
                "operator": "<=",
                "value": 1
            }
        })

    return is_late_stage, findings_json