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
    This function detects a strategy involving silyl protection and deprotection.
    It looks for a silyl group that is present in intermediates but removed in the final product.
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

    has_silyl_intermediate = False
    final_product_has_silyl = False

    def dfs_traverse(node, depth=0):
        nonlocal has_silyl_intermediate, final_product_has_silyl, findings_json

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"]) if node["smiles"] else None

            if mol:
                # Use a robust checker for silyl ether functional group
                # Assuming 'checker' is defined elsewhere or this is a placeholder
                # For the purpose of this refactoring, we'll assume it exists.
                # If 'checker' is not defined, this code will raise a NameError.
                try:
                    has_silyl = checker.has_functional_group(mol, "silyl_ether")
                except NameError:
                    # Placeholder for when 'checker' is not defined in this snippet
                    # In a real scenario, 'checker' would need to be imported or defined.
                    has_silyl = False # Default to False if checker is unavailable

                if has_silyl:
                    findings_json["atomic_checks"]["functional_groups"].append("silyl_ether")

                if depth > 0:  # Intermediate
                    if has_silyl:
                        has_silyl_intermediate = True
                        # Add positional constraint if silyl is found in an intermediate
                        if {"type": "positional", "details": {"target": "silyl_ether", "position": "not_last_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "silyl_ether", "position": "not_last_stage"}})
                else:  # Final product (depth 0)
                    if has_silyl:
                        final_product_has_silyl = True
                    else:
                        # Add negation constraint if silyl is NOT found in the final product
                        if {"type": "negation", "details": {"target": "silyl_ether", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "silyl_ether", "position": "last_stage"}})

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'mol' or 'chemical' type for depth increase
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    result = has_silyl_intermediate and not final_product_has_silyl

    # Return True if silyl groups are present in intermediates but not in the final product
    return result, findings_json
