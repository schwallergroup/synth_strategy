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


# Refactoring for Enumeration: Isolate the list of heterocycles
HETEROCYCLES_OF_INTEREST = ['pyridine', 'triazole', 'piperidine']

def main(route) -> Tuple[bool, Dict]:
    """
    Checks if the final product contains at least two distinct heterocycles
    from the `HETEROCYCLES_OF_INTEREST` list, which includes pyridine, triazole, and piperidine.
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

    final_product_has_multiple_heterocycles = False

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_multiple_heterocycles, findings_json

        # This check should only run on the final product (depth=0)
        if node["type"] == "mol" and depth == 0:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Count how many of the specified heterocycles are present
                found_count = 0
                found_heterocycles = []
                for hetero_name in HETEROCYCLES_OF_INTEREST:
                    # Assuming 'checker' is defined elsewhere or this is a placeholder.
                    # For the purpose of this refactoring, we'll assume it exists.
                    # If 'checker' is not defined, this code would raise an error.
                    # We are only modifying the depth calculation logic.
                    if checker.has_group(mol, hetero_name):
                        found_count += 1
                        found_heterocycles.append(hetero_name)
                
                if found_count >= 2:
                    final_product_has_multiple_heterocycles = True
                    findings_json["atomic_checks"]["ring_systems"].extend(found_heterocycles)
                    # Add the structural constraint when the condition is met
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "distinct_ring_systems_from_list",
                            "target_list": [
                                "pyridine",
                                "triazole",
                                "piperidine"
                            ],
                            "operator": ">=",
                            "value": 2,
                            "position": "last_stage"
                        }
                    })

        # Continue traversing to maintain the original control flow structure,
        # even though the logic only depends on the root node.
        for child in node.get("children", []):
            # Optimization: stop traversal once the condition is met.
            # This is a safe modification as it doesn't change the outcome.
            if final_product_has_multiple_heterocycles:
                break
            
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is not 'reaction' (e.g., 'chemical' or 'mol')
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return final_product_has_multiple_heterocycles, findings_json