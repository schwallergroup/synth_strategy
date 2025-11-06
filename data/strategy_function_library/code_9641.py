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


IMINE_FORMATION_REACTIONS = [
    "Addition of primary amines to aldehydes/thiocarbonyls",
    "Addition of primary amines to ketones/thiocarbonyls",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage imine formation strategy by checking for specific reaction types
    in the final synthetic step. The reactions, defined in the IMINE_FORMATION_REACTIONS
    list, include the addition of primary amines to aldehydes and ketones.
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

    # Track if we found the strategy
    found_imine_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_imine_formation, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Late stage (final step)
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            try:
                # Check if this is one of the target imine formation reactions
                is_imine_formation = False
                for r_name in IMINE_FORMATION_REACTIONS:
                    if checker.check_reaction(r_name, rsmi):
                        is_imine_formation = True
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        break

                if is_imine_formation:
                    found_imine_formation = True
                    # Add the structural constraint if the condition is met
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": [
                                "Addition of primary amines to aldehydes/thiocarbonyls",
                                "Addition of primary amines to ketones/thiocarbonyls"
                            ],
                            "position": "last_stage"
                        }
                    })
            except Exception as e:
                print(f"Error processing reaction SMILES at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            # Depth remains the same when traversing from reaction to chemical
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return found_imine_formation, findings_json
