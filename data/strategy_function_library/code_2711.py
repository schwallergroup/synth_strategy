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
    This function detects a synthetic strategy where a trifluoromethyl group is preserved
    throughout the synthesis.
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

    trifluoromethyl_present_at_all_steps = True
    cf3_pattern = Chem.MolFromSmarts("[C]([F])([F])[F]")

    def dfs_traverse(node, depth=0):
        nonlocal trifluoromethyl_present_at_all_steps, findings_json

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for trifluoromethyl group
                if mol.HasSubstructMatch(cf3_pattern):
                    # If CF3 is present, record it as a functional group found
                    if "trifluoromethyl" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("trifluoromethyl")
                else:
                    # If any non-reagent molecule doesn't have CF3, set flag to False.
                    # This now correctly includes starting materials.
                    if not node.get("in_stock", False):
                        trifluoromethyl_present_at_all_steps = False
                        # Record the structural constraint violation if CF3 is absent where it should be present
                        findings_json["structural_constraints"].append({
                            "type": "negation",
                            "details": {
                                "target": "absence_of_trifluoromethyl_group",
                                "scope": "all_non_stock_molecules"
                            }
                        })

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node.get("type") == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'mol' or 'chemical' type for non-reaction nodes
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return trifluoromethyl_present_at_all_steps, findings_json
