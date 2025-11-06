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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy where a morpholine amide functional group is retained
    throughout multiple steps of the synthesis.
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

    result = False
    morpholine_amide_depths = []

    # SMARTS pattern for morpholine amide
    morpholine_amide_pattern = Chem.MolFromSmarts("[#6](=[O])[#7]1[#6][#6][#8][#6][#6]1")

    def dfs_traverse(node, depth=0):
        nonlocal morpholine_amide_depths, findings_json
        node["depth"] = depth # Assign depth to the current node

        if node["type"] == "mol":
            if node["smiles"]:
                mol = Chem.MolFromSmiles(node["smiles"])

                if mol and mol.HasSubstructMatch(morpholine_amide_pattern):
                    morpholine_amide_depths.append(depth)
                    if "morpholine amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("morpholine amide")

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if morpholine amide is present at multiple unique depths
    if len(set(morpholine_amide_depths)) >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "unique_steps_with_morpholine_amide",
                "operator": ">=",
                "value": 2
            }
        })
        # Sort depths to check range. Sorting the original list with duplicates is fine
        # as it will correctly place the true min and max values at the ends.
        morpholine_amide_depths.sort()

        # Check if the morpholine amide spans a significant portion of the synthesis
        depth_range = morpholine_amide_depths[-1] - morpholine_amide_depths[0]
        if depth_range >= 2:  # Present across at least 3 steps (e.g., depths 1, 2, 3)
            result = True
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "depth_span_of_morpholine_amide",
                    "operator": ">=",
                    "value": 2
                }
            })

    return result, findings_json