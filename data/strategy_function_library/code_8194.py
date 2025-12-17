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


REDUCTIVE_AMINATION_TYPES = [
    "Reductive amination with aldehyde",
    "Reductive amination with ketone",
    "Reductive amination with alcohol",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis employs a late-stage (final or penultimate step) reductive amination.
    The check specifically identifies named reactions for reductive aminations involving aldehydes, ketones, or alcohols,
    as defined in the REDUCTIVE_AMINATION_TYPES list.
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

    found_late_stage_reductive_amination = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_reductive_amination, findings_json

        # For reaction nodes, check if it's a late-stage reductive amination
        if node["type"] == "reaction" and depth <= 2:  # Late stage (low depth)
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Direct check for reductive amination reactions
                for rxn_type in REDUCTIVE_AMINATION_TYPES:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Found late-stage reductive amination at depth {depth}: {rsmi}")
                        found_late_stage_reductive_amination = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        
                        # Add structural constraint if not already added
                        positional_constraint = {
                            "type": "positional",
                            "details": {
                                "targets": [
                                    "Reductive amination with aldehyde",
                                    "Reductive amination with ketone",
                                    "Reductive amination with alcohol"
                                ],
                                "position_description": "Occurs within the last three reaction steps (depth <= 2)."
                            }
                        }
                        if positional_constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(positional_constraint)
                        
                        # No need to break here, as we want to record all matching types if multiple apply

        # Process children with incremented depth
        for child in node.get("children", []):
            # New logic: depth increases only from chemical to reaction node
            if node["type"] == "chemical":
                dfs_traverse(child, depth + 1)
            else: # node["type"] == "reaction"
                dfs_traverse(child, depth)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage reductive amination detected: {found_late_stage_reductive_amination}")
    return found_late_stage_reductive_amination, findings_json
