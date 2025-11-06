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
    This function detects a strategy involving the use of a difluoro acetal protecting group
    that is maintained throughout the synthesis.
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

    depths_with_difluoro_acetal = set()
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal depths_with_difluoro_acetal, max_depth, findings_json

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[O][C]([F])([F])[O]")):
                    depths_with_difluoro_acetal.add(depth)
                    if "difluoro acetal" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("difluoro acetal")
            except Exception as e:
                pass

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if difluoro acetal is present at multiple depths (maintained throughout synthesis)
    result = len(depths_with_difluoro_acetal) >= 3 and max(depths_with_difluoro_acetal) >= 3

    if len(depths_with_difluoro_acetal) >= 3:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "stages_with_difluoro_acetal",
                "operator": ">=",
                "value": 3
            }
        })

    return result, findings_json
