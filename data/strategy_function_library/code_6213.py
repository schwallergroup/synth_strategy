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
    Detects if a trifluoromethyl-substituted pyridine scaffold is maintained throughout the synthesis.
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

    # Track molecules at each depth
    molecules_by_depth = {}
    max_depth = -1
    
    result = True # Initialize the overall result flag

    def dfs_traverse(node, depth=0):
        nonlocal max_depth

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "mol" and "smiles" in node:
            if depth not in molecules_by_depth:
                molecules_by_depth[depth] = []
            molecules_by_depth[depth].append(node["smiles"])

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not 'reaction' (e.g., 'chemical' or 'mol')
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if trifluoromethyl-pyridine is present at all depths
    for depth in range(max_depth + 1):
        if depth not in molecules_by_depth:
            result = False
            break

        found_at_depth = False
        for smiles in molecules_by_depth[depth]:
            # A single substructure search is more robust and efficient
            # than checking for groups separately and then verifying connectivity.
            # SMARTS: aromatic pyridine ('c1ncccc1') attached to a CF3 group ('C(F)(F)F').
            if checker.check_substructure('c1(ncccc1)C(F)(F)F', smiles):
                found_at_depth = True
                # Record atomic checks when the combined substructure is found
                if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("pyridine")
                if "trifluoromethyl" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("trifluoromethyl")
                break

        if not found_at_depth:
            result = False
            break
    
    # If the overall result is True, add the structural constraint
    if result:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "trifluoromethyl-pyridine",
                "position": "all_stages"
            }
        })

    return result, findings_json
