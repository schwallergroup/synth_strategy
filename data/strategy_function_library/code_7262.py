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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a linear synthesis strategy where an indole scaffold is
    present in molecules across multiple steps. It invalidates the strategy if a
    convergent step involving multiple indole-containing reactants is found.
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

    # Initialize tracking variables
    reaction_count = 0
    molecules_with_methoxy_indole = []
    result = True # Overall result flag

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, result, findings_json, molecules_with_methoxy_indole

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check for indole ring
            has_indole = checker.check_ring("indole", mol_smiles)

            if has_indole:
                findings_json["atomic_checks"]["ring_systems"].append("indole")
                # Store molecules with indole along with their depth
                molecules_with_methoxy_indole.append((mol_smiles, depth))

        elif node["type"] == "reaction":
            reaction_count += 1

            # Check if this is a convergent step (more than one reactant with indole)
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                # Count reactants containing indole
                indole_reactants = [r for r in reactants if checker.check_ring("indole", r)]
                if len(indole_reactants) > 1:
                    result = False # Set overall result to False
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "convergent_indole_reaction",
                            "description": "A reaction step must not have more than one reactant containing an indole ring."
                        }
                    })
                    return False

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            traverse_result = dfs_traverse(child, next_depth)
            if traverse_result is False:  # If any child returns False, propagate it up
                result = False # Propagate overall result
                return False

        return True

    # Start traversal from the root
    linear_synthesis_traversal_result = dfs_traverse(route)

    # Sort molecules by depth (ascending order)
    molecules_with_methoxy_indole.sort(key=lambda x: x[1])

    # Check if we have at least 3 reaction steps
    has_sufficient_steps = reaction_count >= 3
    if has_sufficient_steps:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction",
                "operator": ">=",
                "value": 3,
                "description": "The total number of reaction steps in the route must be at least 3."
            }
        })

    # Check if we have indole in at least 3 molecules (start, intermediate, end)
    has_methoxy_indole_throughout = len(molecules_with_methoxy_indole) >= 3
    if has_methoxy_indole_throughout:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "indole_presence",
                "operator": ">=",
                "value": 3,
                "description": "The number of molecules containing an indole ring throughout the route must be at least 3."
            }
        })

    # Check if the strategy is present
    strategy_present = linear_synthesis_traversal_result and has_sufficient_steps and has_methoxy_indole_throughout

    # Update the overall result based on all conditions
    result = result and strategy_present

    print(
        f"Strategy detection results: methoxy_indole_count={len(molecules_with_methoxy_indole)}, "
        f"linear_synthesis={linear_synthesis_traversal_result}, reaction_count={reaction_count}"
    )

    return result, findings_json
