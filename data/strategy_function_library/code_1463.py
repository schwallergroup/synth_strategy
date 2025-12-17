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
    Detects if an indole heterocycle and a fluorine atom each appear in at least two separate molecules within the synthesis tree.
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

    indole_depths = []
    fluorine_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal indole_depths, fluorine_depths, findings_json

        if node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]

            # Check for indole core using the checker function
            if checker.check_ring("indole", mol_smiles):
                indole_depths.append(depth)
                if "indole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("indole")

            # Check for fluorine atoms using the checker function
            if checker.has_atom("F", mol_smiles):
                fluorine_depths.append(depth)
                if "fluorine" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("fluorine")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # Depth increases only if not a reaction node
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if both motifs are present at multiple depths (preserved throughout)
    indole_preserved = len(indole_depths) >= 2
    fluorine_preserved = len(fluorine_depths) >= 2
    both_preserved = indole_preserved and fluorine_preserved

    if indole_preserved:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "indole",
                "operator": ">=",
                "value": 2
            }
        })
    if fluorine_preserved:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "fluorine",
                "operator": ">=",
                "value": 2
            }
        })
    if both_preserved:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "indole",
                    "fluorine"
                ]
            }
        })

    return both_preserved, findings_json
