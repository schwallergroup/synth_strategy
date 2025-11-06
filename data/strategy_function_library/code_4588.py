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


def main(route) -> Tuple[bool, Dict]:
    """This function detects the reduction of a nitro group to a primary amine, a common tactic for using the nitro group as a masked amine."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    has_nitro_reduction = False

    # The inner function signature is updated to accept standard context.
    def dfs_traverse(node, reaction, depth, max_depth):
        nonlocal has_nitro_reduction, findings_json

        if reaction:
            # The original implementation was flawed, checking only for the presence/absence
            # of functional groups, which leads to false positives. It is replaced by a
            # single, robust checker that verifies the specific functional group conversion.
            if checker.is_fg_conversion(reaction, "nitro", "primary_amine"):
                has_nitro_reduction = True
                findings_json["atomic_checks"]["named_reactions"].append("nitro_reduction")
                findings_json["atomic_checks"]["functional_groups"].append("nitro")
                findings_json["atomic_checks"]["functional_groups"].append("primary_amine")

        # The wrapper's traversal logic is updated to pass context to children.
        for child in node.get("children", []):
            child_reaction = checker.get_reaction(child)
            
            # Determine the new depth based on the current node's type
            # Assuming node.get('type') returns 'chemical' or 'reaction'
            current_node_type = node.get('type')
            if current_node_type == 'reaction':
                # Depth remains the same when traversing from a reaction node to a chemical node
                new_depth = depth
            else: # Assuming 'chemical' or any other type that implies a chemical node
                # Depth increases when traversing from a chemical node to a reaction node
                new_depth = depth + 1

            dfs_traverse(child, child_reaction, new_depth, max_depth)

    # The initial call is updated to provide the starting context.
    # We assume helper functions exist to get the initial reaction and max_depth.
    initial_reaction = checker.get_reaction(route)
    max_depth = checker.get_max_depth(route)
    dfs_traverse(route, initial_reaction, 1, max_depth)
    
    if has_nitro_reduction:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "nitro_reduction",
                "operator": ">=",
                "value": 1
            }
        })

    return has_nitro_reduction, findings_json