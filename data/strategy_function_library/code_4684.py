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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves a late-stage nitrile hydrolysis (C≡N → CONH₂)
    in the final step of the synthesis.
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

    nitrile_hydrolysis_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_hydrolysis_found, findings_json

        # For molecule nodes, check if it's the target molecule (depth 0)
        if node["type"] == "mol" and depth == 0:
            # Check if any of the child reactions is a nitrile hydrolysis
            for child in node.get("children", []):
                if nitrile_hydrolysis_found:
                    break
                if child["type"] == "reaction" and child.get("metadata", {}).get("rsmi"):
                    rsmi = child["metadata"]["rsmi"]
                    if checker.check_reaction("Nitrile to amide", rsmi):
                        nitrile_hydrolysis_found = True
                        findings_json["atomic_checks"]["named_reactions"].append("Nitrile to amide")
                        # Add the structural constraint if the condition is met
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "Nitrile to amide",
                                "position": "last_stage"
                            }
                        })

        # Traverse children with incremented depth
        if not nitrile_hydrolysis_found:
            for child in node.get("children", []):
                # New logic for depth calculation
                new_depth = depth
                if node["type"] != "reaction":  # If current node is chemical, depth increases
                    new_depth = depth + 1
                # If current node is reaction, depth remains the same
                dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return nitrile_hydrolysis_found, findings_json