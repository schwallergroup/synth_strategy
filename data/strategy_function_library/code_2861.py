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
    Detects if the synthesis uses Boc protection strategy throughout most steps
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

    boc_protected_intermediates = 0
    total_intermediates = 0
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_protected_intermediates, total_intermediates, findings_json

        if node["type"] == "mol" and "smiles" in node and depth > 0:  # Skip final product
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                total_intermediates += 1

                # Check for Boc group using the checker function
                # Assuming 'checker' is available in the context or passed as an argument
                # For this refactoring, we assume 'checker' is globally accessible or mocked for testing.
                # In a real scenario, 'checker' would need to be defined or passed.
                try:
                    if checker.has_functional_group(mol, "boc"):
                        boc_protected_intermediates += 1
                        if "boc" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("boc")
                except NameError:
                    # Handle case where 'checker' is not defined for this snippet
                    pass # Or log an error

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # Assuming 'mol' or other types that should increase depth
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if more than half of intermediates are Boc-protected
    if total_intermediates > 0:
        ratio = boc_protected_intermediates / total_intermediates
        if ratio >= 0.5:
            result = True
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "proportion_of_intermediates_with_boc",
                    "operator": ">=",
                    "value": 0.5
                }
            })

    return result, findings_json