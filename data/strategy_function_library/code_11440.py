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
    Checks if any of the starting materials (leaf nodes in the synthesis tree) is a 
    poly-fluorinated aromatic compound, defined as having two or more fluorine atoms 
    attached to an aromatic ring. This strategy identifies routes that begin with 
    these key building blocks.
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

    fluorinated_aromatic_found = False

    # Helper to determine the maximum depth of the synthesis tree.
    def _get_max_depth(node, depth=1):
        if not node.get("children", []):
            return depth
        return max(_get_max_depth(child, depth + 1) for child in node["children"])

    max_depth = _get_max_depth(route) if route else 0

    def dfs_traverse(node, depth):
        nonlocal fluorinated_aromatic_found, findings_json
        if fluorinated_aromatic_found:
            return

        # A "building block" is a starting material, which is a leaf node.
        # In a depth-first traversal from the root, leaf nodes are at max_depth.
        if node["type"] == "mol" and "smiles" in node and depth == max_depth:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol is not None:
                # Check for at least two fluorines on an aromatic ring
                matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[c][F]"))
                if len(matches) >= 2:
                    fluorinated_aromatic_found = True
                    findings_json["atomic_checks"]["functional_groups"].append("aryl fluoride")
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "aryl fluoride",
                            "operator": ">=",
                            "value": 2
                        }
                    })
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "aryl fluoride",
                            "position": "first_stage"
                        }
                    })
                    return

        # Traverse children, incrementing depth
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'mol' (chemical)
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root (depth=1)
    if route:
        dfs_traverse(route, 1)
    return fluorinated_aromatic_found, findings_json