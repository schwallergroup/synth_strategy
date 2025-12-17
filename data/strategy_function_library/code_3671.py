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
    Detects if the synthesis follows a primarily linear strategy with a late-stage convergent step.
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

    # Track the number of linear vs. convergent steps
    step_count = 0
    convergent_steps = 0
    late_convergent_step = False

    def dfs_traverse(node, depth=0):
        nonlocal step_count, convergent_steps, late_convergent_step, findings_json

        if node["type"] == "reaction":
            step_count += 1
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                # If multiple reactants, it's a convergent step
                if len(reactants) > 1:
                    convergent_steps += 1
                    findings_json["atomic_checks"]["named_reactions"].append("convergent_step")
                    # If depth is 0 or 1, it's a late-stage convergent step
                    if depth <= 1:
                        late_convergent_step = True

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth when going from chemical to reaction
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # If most steps are linear but there's at least one late convergent step
    result = late_convergent_step and (convergent_steps <= step_count / 2)

    if late_convergent_step:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "convergent_step",
                "position": "last_two_stages"
            }
        })
    
    if convergent_steps <= step_count / 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "convergent_step",
                "operator": "<=",
                "value": "total_reaction_steps / 2"
            }
        })

    return result, findings_json