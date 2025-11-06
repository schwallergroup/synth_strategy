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


SONOGASHIRA_REACTION_TYPES = [
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl OTf",
    "Sonogashira alkyne_aryl OTf",
    "Sonogashira acetylene_alkenyl halide",
    "Sonogashira alkyne_alkenyl halide",
    "Sonogashira acetylene_alkenyl OTf",
    "Sonogashira alkyne_alkenyl OTf",
    "Sonogashira acetylene_acyl halide",
    "Sonogashira alkyne_acyl halide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for late-stage Sonogashira couplings using a predefined list of reaction types. This includes couplings where the electrophile is an aryl, alkenyl, or acyl halide/triflate.
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

    sonogashira_found = False

    def dfs_traverse(node, depth=0):
        nonlocal sonogashira_found, findings_json

        if node["type"] == "reaction" and depth <= 2:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for Sonogashira coupling using the checker function
                for reaction_type in SONOGASHIRA_REACTION_TYPES:
                    if checker.check_reaction(reaction_type, rsmi):
                        sonogashira_found = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # Only increase depth if current node is not a reaction
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    if sonogashira_found:
        # Add the structural constraint if a Sonogashira reaction was found within the depth limit
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Sonogashira coupling",
                "position": "within_last_3_stages"
            }
        })

    return sonogashira_found, findings_json
