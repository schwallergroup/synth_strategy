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
    This function detects a synthetic strategy involving a multi-component reaction
    in the early stage of synthesis (combining 3+ fragments).
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

    found_multi_component_reaction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_multi_component_reaction, findings_json

        node["depth"] = depth # Assign depth to the current node

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Count valid reactants (non-empty SMILES)
            valid_reactants = [smi for smi in reactants_smiles if smi]

            # Check if this is a multi-component reaction (3+ reactants)
            if len(valid_reactants) >= 3:
                findings_json["structural_constraints"].append({
                    "type": "count",
                    "details": {
                        "target": "reactants_in_reaction",
                        "operator": ">=",
                        "value": 3
                    }
                })
                # Check if this is not the final step (depth 2+)
                if node.get("depth", 0) >= 2:
                    found_multi_component_reaction = True
                    findings_json["atomic_checks"]["named_reactions"].append("multi-component reaction")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "multi-component reaction",
                            "position": "not_final_or_penultimate_stage"
                        }
                    })

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_multi_component_reaction, findings_json
