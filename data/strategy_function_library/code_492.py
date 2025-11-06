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

HUISGEN_REACTION_VARIANTS = [
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen_Cu-catalyzed_1,4-subst",
    "Huisgen_Ru-catalyzed_1,5_subst",
    "Huisgen 1,3 dipolar cycloaddition",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route uses a late-stage Huisgen cycloaddition by checking for specific named reaction variants. The checked variants include the thermal, copper-catalyzed (CuAAC), and ruthenium-catalyzed (RuAAC) versions.
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

    found_click_chemistry = False

    def dfs_traverse(node, depth=0):
        nonlocal found_click_chemistry, findings_json

        if node["type"] == "reaction":
            try:
                # Check if this is a late-stage reaction (depth <= 1)
                is_late_stage = depth <= 1

                if is_late_stage:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    
                    # Check if this is a Huisgen cycloaddition reaction
                    is_huisgen_reaction = False
                    for name in HUISGEN_REACTION_VARIANTS:
                        if checker.check_reaction(name, rsmi):
                            is_huisgen_reaction = True
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                            break

                    if is_huisgen_reaction:
                        found_click_chemistry = True
                        # Add the structural constraint if a late-stage Huisgen reaction is found
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "targets": [
                                    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                                    "Huisgen_Cu-catalyzed_1,4-subst",
                                    "Huisgen_Ru-catalyzed_1,5_subst",
                                    "Huisgen 1,3 dipolar cycloaddition"
                                ],
                                "position": "late_stage",
                                "condition": "depth <= 1"
                            }
                        })
            except Exception:
                # Errors in reaction analysis are silently ignored to avoid crashing the traversal.
                pass

        # Traverse all children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return found_click_chemistry, findings_json
