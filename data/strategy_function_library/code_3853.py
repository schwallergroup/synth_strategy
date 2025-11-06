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
    Detects if a fluorinated aryl group is maintained throughout the synthesis.
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

    fluoro_aryl_depths = set()

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Fluorinated aryl SMARTS pattern
                    fluoro_aryl_pattern = Chem.MolFromSmarts("[c][F]")
                    if mol.HasSubstructMatch(fluoro_aryl_pattern):
                        fluoro_aryl_depths.add(depth)
                        if "fluorinated aryl" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("fluorinated aryl")
            except:
                pass

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'mol' or 'chemical' type for other nodes
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if fluorinated aryl is present at multiple depths
    result = len(fluoro_aryl_depths) >= 2

    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "presence_of_fluorinated_aryl_at_distinct_depths",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json