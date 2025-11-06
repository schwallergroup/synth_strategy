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


import re
from rdkit import Chem

# Assume 'checker' is a pre-existing library with the following functions:
# checker.is_fg_formed(reaction_object, fg_name) -> bool

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a trifluoromethyl group is introduced in a late-stage step of the synthesis, defined as occurring within the final two steps (depth 1 or 2).
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

    cf3_introduced = False
    introduction_depth = -1

    def dfs_traverse(node, depth):
        nonlocal cf3_introduced, introduction_depth, findings_json

        # Process the node if it's a reaction and has a reaction object.
        if node.get("type") == "reaction" and "reaction" in node:
            # The check for CF3 formation is now a single, robust call,
            # replacing manual, error-prone SMARTS/SMILES parsing.
            if checker.is_fg_formed(node["reaction"], "trifluoromethyl"):
                cf3_introduced = True
                findings_json["atomic_checks"]["functional_groups"].append("trifluoromethyl")
                findings_json["atomic_checks"]["named_reactions"].append("trifluoromethyl_formation") # Conceptual reaction name
                # This logic correctly finds the latest-stage (smallest depth)
                # introduction by overwriting depth on each find as the DFS unwinds
                # from its deepest point.
                introduction_depth = depth

        # Recursively traverse to children (reactants), incrementing the depth.
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction.
            # Depth remains the same when going from reaction to chemical.
            if node.get("type") == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' type or other non-reaction type
                dfs_traverse(child, depth + 1)

    # Start the depth-first search from the root of the synthesis tree (the final product).
    # The final step is at depth=1.
    dfs_traverse(route, 1)

    # The strategy is positive if a CF3 group was introduced at a late stage (depth 1 or 2).
    result = cf3_introduced and introduction_depth <= 2

    if result:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "trifluoromethyl_formation",
                "position": "last_two_stages"
            }
        })

    return result, findings_json