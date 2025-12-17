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


AZIDE_ALKYNE_CYCLOADDITIONS = [
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen_Cu-catalyzed_1,4-subst",
    "Huisgen_Ru-catalyzed_1,5_subst",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves one of the specified azide-alkyne cycloaddition reactions
    (e.g., Huisgen thermal, Cu-catalyzed, or Ru-catalyzed).
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

    click_chemistry_used = False

    def dfs_traverse(node):
        nonlocal click_chemistry_used, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for any of the defined azide-alkyne cycloaddition reactions
            for reaction_name in AZIDE_ALKYNE_CYCLOADDITIONS:
                if checker.check_reaction(reaction_name, rsmi):
                    print(f"Click chemistry detected: {reaction_name}")
                    click_chemistry_used = True
                    if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    break # Found one, no need to check others for this node

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if click_chemistry_used:
        # Add the structural constraint if any of the target reactions were found
        findings_json["structural_constraints"].append(
            {
                "type": "count",
                "details": {
                    "target": [
                        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                        "Huisgen_Cu-catalyzed_1,4-subst",
                        "Huisgen_Ru-catalyzed_1,5_subst"
                    ],
                    "operator": ">=",
                    "value": 1
                }
            }
        )

    return click_chemistry_used, findings_json
