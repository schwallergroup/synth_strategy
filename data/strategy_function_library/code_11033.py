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
    This function detects a synthetic strategy involving a late-stage Sonogashira coupling
    (terminal alkyne + aryl halide) in the final step of a synthesis.
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

    sonogashira_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal sonogashira_detected, findings_json

        if node["type"] == "reaction" and depth == 1:  # Final reaction step
            try:
                if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]

                    # Check if this is a Sonogashira coupling reaction directly
                    for rxn_type in SONOGASHIRA_REACTION_TYPES:
                        if checker.check_reaction(rxn_type, rsmi):
                            sonogashira_detected = True
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                            # Add the structural constraint if detected in the final step
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "targets": [
                                        "Sonogashira acetylene_aryl halide",
                                        "Sonogashira alkyne_aryl halide",
                                        "Sonogashira acetylene_aryl OTf",
                                        "Sonogashira alkyne_aryl OTf",
                                        "Sonogashira acetylene_alkenyl halide",
                                        "Sonogashira alkyne_alkenyl halide",
                                        "Sonogashira acetylene_alkenyl OTf",
                                        "Sonogashira alkyne_alkenyl OTf",
                                        "Sonogashira acetylene_acyl halide",
                                        "Sonogashira alkyne_acyl halide"
                                    ],
                                    "position": "last_stage"
                                }
                            })
                            return

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversing
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

    return sonogashira_detected, findings_json
