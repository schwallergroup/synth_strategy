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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the final step in a synthesis involves
    oxidation of a sulfide to a sulfoxide.
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

    found_late_stage_oxidation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_oxidation, findings_json

        # The final product is at depth 0 (mol node)
        # The final reaction step is at depth 1 (reaction node)
        if (
            depth == 1
            and node["type"] == "reaction"
            and "metadata" in node
            and "mapped_reaction_smiles" in node["metadata"]
        ):
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            try:
                # Check if this is a sulfide to sulfoxide oxidation
                if checker.check_reaction("Sulfanyl to sulfinyl", rsmi):
                    print(f"Found late-stage sulfide oxidation via reaction check: {rsmi}")
                    found_late_stage_oxidation = True
                    findings_json["atomic_checks"]["named_reactions"].append("Sulfanyl to sulfinyl")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Sulfanyl to sulfinyl",
                            "position": "last_stage"
                        }
                    })
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage sulfoxide formation detection result: {found_late_stage_oxidation}")
    return found_late_stage_oxidation, findings_json
