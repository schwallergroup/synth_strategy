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


# Refactoring for Enumeration: The list of reaction types is moved to the module level.
NITROGEN_ATTACHMENT_REACTION_TYPES = [
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Buchwald-Hartwig",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Goldberg coupling",
    "Goldberg coupling aryl amine-aryl chloride",
    "Goldberg coupling aryl amide-aryl chloride",
    "Reductive amination",
    "Reductive amination with aldehyde",
    "Reductive amination with ketone",
    "Reductive amination with alcohol",
    "Amination",
    "Ullmann-Goldberg Substitution amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy involving three or more sequential nitrogen attachment steps.
    A step is identified as a nitrogen attachment if it matches a reaction type
    from the predefined list NITROGEN_ATTACHMENT_REACTION_TYPES.
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

    n_attachment_steps = 0
    sequential = True
    last_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal n_attachment_steps, sequential, last_depth, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            try:
                is_n_attachment_reaction = False
                for rxn_type in NITROGEN_ATTACHMENT_REACTION_TYPES:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_n_attachment_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if is_n_attachment_reaction:
                    n_attachment_steps += 1

                    # Check if steps are sequential in the synthetic direction
                    if last_depth != -1:
                        if depth <= last_depth:  # Non-sequential in retrosynthetic direction
                            sequential = False
                    last_depth = depth
            except Exception:
                # Silently ignore checker errors for robustness
                pass

        # Process children only if the sequence is still valid
        if sequential:
            for child in node.get("children", []):
                # New logic for depth calculation
                new_depth = depth
                if node["type"] != "reaction":  # If current node is chemical, depth increases
                    new_depth = depth + 1
                dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = n_attachment_steps >= 3 and sequential

    # Add structural constraints if the conditions are met
    if n_attachment_steps >= 3:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "nitrogen_attachment_reaction",
                "operator": ">=",
                "value": 3
            }
        })
    if sequential:
        # This constraint is met if 'sequential' remains True throughout the traversal
        # and at least two attachment steps were found to establish a sequence.
        # The original JSON implies a sequence of at least two, which is implicitly handled
        # by 'sequential' remaining True if multiple steps are found.
        # We add it if the overall 'sequential' condition holds.
        if n_attachment_steps >= 2: # Need at least two steps to form a sequence
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "targets": [
                        "nitrogen_attachment_reaction",
                        "nitrogen_attachment_reaction"
                    ]
                }
            })

    return result, findings_json
