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
    Detects if the synthesis follows a linear strategy with a late-stage coupling of fragments.
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

    late_coupling = False
    linear_steps = 0
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal late_coupling, linear_steps, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if this is a coupling reaction (combining multiple fragments)
                if len(reactants) >= 2 and depth <= 1:  # Late stage
                    late_coupling = True
                    # Record the atomic check for fragment_coupling
                    if "fragment_coupling" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("fragment_coupling")

            linear_steps += 1

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # A synthesis is considered linear if it has multiple steps and ends with a coupling
    result = late_coupling and linear_steps >= 3 and max_depth >= 3

    # Record structural constraints if conditions are met
    if late_coupling:
        # This corresponds to the positional constraint: "target": "fragment_coupling", "position": "in_final_two_steps"
        # The 'depth <= 1' check in dfs_traverse already implies 'late-stage' or 'in_final_two_steps'
        if {"type": "positional", "details": {"target": "fragment_coupling", "position": "in_final_two_steps"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "fragment_coupling", "position": "in_final_two_steps"}})

    if linear_steps >= 3:
        # This corresponds to the count constraint: "target": "any_reaction", "operator": ">=", "value": 3
        if {"type": "count", "details": {"target": "any_reaction", "operator": ">=", "value": 3}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "any_reaction", "operator": ">=", "value": 3}})

    if max_depth >= 3:
        # This corresponds to the count constraint: "target": "synthesis_depth", "operator": ">=", "value": 3
        if {"type": "count", "details": {"target": "synthesis_depth", "operator": ">=", "value": 3}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "synthesis_depth", "operator": ">=", "value": 3}})

    return result, findings_json
