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


# Refactoring for Enumeration: Isolate the list of reaction types
CH_FUNCTIONALIZATION_REACTIONS = [
    "Aromatic chlorination",
    "Aromatic bromination",
    "Aromatic iodination",
    "Aromatic fluorination",
    "Friedel-Crafts acylation",
    "Directed ortho metalation of arenes",
    "Minisci (para)",
    "Minisci (ortho)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects an early-stage C-H functionalization strategy. This is defined as a reaction
    occurring within the first three steps of the synthesis (i.e., depth >= max_depth - 2)
    that matches a specific, known C-H functionalization reaction type.
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

    has_early_ch_functionalization = False

    # Helper to find max_depth, necessary for robust 'early-stage' definition
    def _get_max_depth(node, current_depth=0):
        if not node.get("children"):
            return current_depth
        # The depth should only increase when going from a chemical to a reaction node.
        # When going from a reaction to a chemical node, the depth remains the same.
        next_depth = current_depth + 1 if node.get("type") != "reaction" else current_depth
        return max(_get_max_depth(child, next_depth) for child in node["children"])

    max_depth = _get_max_depth(route)

    # Propagating Context: The signature is updated to include max_depth
    def dfs_traverse(node, depth, max_depth):
        nonlocal has_early_ch_functionalization, findings_json

        # Conditional Modification: Use max_depth for a robust 'early-stage' check
        if node.get("type") == "reaction" and depth >= max_depth - 2:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                
                # Code Removal: The flawed fallback logic has been removed.
                # The check is now solely based on a robust list of named reactions.
                for reaction_type in CH_FUNCTIONALIZATION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        has_early_ch_functionalization = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        # Add the structural constraint if it's an early-stage reaction
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "targets": [
                                    "Aromatic chlorination",
                                    "Aromatic bromination",
                                    "Aromatic iodination",
                                    "Aromatic fluorination",
                                    "Friedel-Crafts acylation",
                                    "Directed ortho metalation of arenes",
                                    "Minisci (para)",
                                    "Minisci (ortho)"
                                ],
                                "position": "first_three_stages"
                            }
                        })
                        break
            except Exception:
                # Silently ignore errors in malformed nodes
                pass

        # Traverse children, propagating context
        for child in node.get("children", []):
            # New Logic: Depth increases only when traversing from chemical to reaction.
            # Depth remains the same when traversing from reaction to chemical.
            next_depth = depth + 1 if node.get("type") != "reaction" else depth
            dfs_traverse(child, next_depth, max_depth)

    # The initial call is updated to pass the necessary context.
    dfs_traverse(route, 0, max_depth)

    return has_early_ch_functionalization, findings_json
