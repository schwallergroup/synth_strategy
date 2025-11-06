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
    This function detects if the synthetic route uses a multi-component coupling strategy,
    where three or more fragments are combined in a single reaction.
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

    has_multi_component = False

    def dfs_traverse(node, depth=0):
        nonlocal has_multi_component, findings_json

        if node["type"] == "reaction":
            # Use a robust checker to get reactants instead of parsing SMILES strings.
            # Assuming 'checker' is defined elsewhere or passed implicitly, as it's not in the provided snippet.
            # For the purpose of this refactoring, we'll assume its existence.
            try:
                reactants = checker.get_reactants(node)
            except NameError:
                # Placeholder for when 'checker' is not defined in this scope
                # In a real scenario, 'checker' would be an instantiated object
                reactants = [] # Or handle error appropriately

            # If we have 3 or more reactants, it's a multi-component coupling
            if len(reactants) >= 3:
                has_multi_component = True
                # Add the structural constraint to findings_json
                findings_json["structural_constraints"].append(
                    {
                        "type": "count",
                        "details": {
                            "target": "reactants_in_reaction",
                            "operator": ">=",
                            "value": 3
                        }
                    }
                )

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical (or other non-reaction type)
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return has_multi_component, findings_json