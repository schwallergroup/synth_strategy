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

def main(route) -> Tuple[bool, Dict]:
    """This function detects if a fluorination reaction, including aromatic fluorination, is used at any step in the synthesis."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    early_fluorination = False

    def dfs_traverse(node, depth=0):
        nonlocal early_fluorination, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Rely on robust checkers to identify a fluorination reaction.
                # This is more reliable than manual checks and avoids common pitfalls.
                is_fluorination = checker.check_reaction("Fluorination", rsmi)
                is_aromatic_fluorination = checker.check_reaction("Aromatic fluorination", rsmi)

                if is_fluorination:
                    early_fluorination = True
                    if "Fluorination" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Fluorination")
                
                if is_aromatic_fluorination:
                    early_fluorination = True
                    if "Aromatic fluorination" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Aromatic fluorination")

            except (KeyError, IndexError):
                # Handle cases where rsmi is missing or malformed.
                pass

        # Continue traversal to check all reactions in the route.
        # An early exit could be more efficient but would alter the original control flow.
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' type or any other type that should increase depth
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if early_fluorination:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": [
                    "Fluorination",
                    "Aromatic fluorination"
                ],
                "operator": ">=",
                "value": 1
            }
        })

    return early_fluorination, findings_json
