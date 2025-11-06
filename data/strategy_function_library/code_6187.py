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
    Detects syntheses where at least one path from a starting material to the final product consists of a convergent final step followed by a linear sequence of reactions. A convergent step is defined as having three or more reactant molecules in the reaction SMILES. A linear step is defined as having one or two reactants. The strategy requires at least three sequential linear steps preceding the final convergent step.

    Args:
        route: A synthesis route JSON object

    Returns:
        bool: True if the route follows a linear-to-convergent strategy
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

    # Track all paths through the synthesis
    synthesis_paths = []
    result = False

    def dfs_traverse(node, path, depth):
        nonlocal findings_json
        if node["type"] == "mol":
            # If this is a leaf node (starting material), we've completed a path
            if node.get("in_stock", False) or not node.get("children", []):
                # Save the completed path
                synthesis_paths.append(path.copy())
                return

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            num_reactants = len(reactants)

            # Add this reaction to the current path with its properties
            path_entry = {
                "depth": depth,
                "num_reactants": num_reactants,
                "is_convergent": num_reactants >= 3,
                "rsmi": rsmi,
            }
            path.append(path_entry)

            if num_reactants >= 3:
                if "convergent_reaction" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("convergent_reaction")
            elif 1 <= num_reactants <= 2:
                if "linear_reaction" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("linear_reaction")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is 'chemical' (or 'mol'), depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, path.copy(), next_depth)

    # Start traversal from the root with initial depth 0
    dfs_traverse(route, [], 0)

    # Analyze each synthesis path for linear-to-convergent pattern
    for path in synthesis_paths:
        # In retrosynthetic analysis, the first reaction in the path is the final step in forward synthesis
        if not path or not path[0]["is_convergent"]:
            continue

        # Count sequential linear steps starting from the second reaction (after the convergent step)
        linear_count = 0
        for i in range(1, len(path)):
            if 1 <= path[i]["num_reactants"] <= 2:
                linear_count += 1
            else:
                # Stop counting at first non-linear step
                break

        # If we found a valid linear-to-convergent path, set result to True and record structural constraint
        if linear_count >= 3:
            result = True
            # Add the structural constraint as defined in the strategy JSON
            findings_json["structural_constraints"].append(
                {
                    "type": "sequence",
                    "details": {
                        "description": "A final convergent step must be immediately preceded by a block of at least 3 sequential linear steps.",
                        "ordered_events": [
                            {
                                "target": "linear_reaction",
                                "constraint": {
                                    "type": "count",
                                    "details": {
                                        "operator": ">=",
                                        "value": 3,
                                        "is_sequential": True
                                    }
                                }
                            },
                            {
                                "target": "convergent_reaction",
                                "constraint": {
                                    "type": "positional",
                                    "details": {
                                        "position": "last_stage"
                                    }
                                }
                            }
                        ]
                    }
                }
            )
            # Since we only need to find one such path, we can break here
            break

    # No valid linear-to-convergent path found
    return result, findings_json
