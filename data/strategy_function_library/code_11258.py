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


COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki coupling",
    "Negishi coupling",
    "Stille reaction",
    "Heck reaction",
    "Sonogashira",
    "Buchwald-Hartwig",
    "N-arylation",
    "Ullmann condensation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis uses a late-stage convergent coupling strategy
    where three or more fragments are combined in the final step.
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

    convergent_coupling_found = False
    final_step_node = None

    # First pass to find the final step (depth 0)
    def find_final_step(node, depth=0):
        nonlocal final_step_node

        # If this is a molecule node and it's the target molecule (depth 0)
        if node["type"] == "mol" and depth == 0:
            # Look for its parent reaction
            for child in node.get("children", []):
                if child["type"] == "reaction":
                    final_step_node = child
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "coupling_reaction",
                            "position": "last_stage"
                        }
                    })
                    break

        # Continue traversing
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            find_final_step(child, new_depth)

    # Function to check if a reaction is a coupling reaction
    def is_coupling_reaction(rsmi):
        # Check common coupling reactions
        for rxn_type in COUPLING_REACTIONS_OF_INTEREST:
            if checker.check_reaction(rxn_type, rsmi):
                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                return True

        return False

    # Second pass to check for convergent coupling
    def check_convergent_coupling(node):
        nonlocal convergent_coupling_found

        if node == final_step_node:
            rsmi = node["metadata"].get("rsmi", "")
            if rsmi:
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                # If there are 3 or more reactants in the final step, check if it's a coupling
                if len(reactants) >= 3:
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "reactants_in_final_reaction",
                            "operator": ">=",
                            "value": 3
                        }
                    })
                    if is_coupling_reaction(rsmi):
                        convergent_coupling_found = True

        # Continue traversing
        for child in node.get("children", []):
            check_convergent_coupling(child)

    # Start traversal from the root to find final step
    find_final_step(route)

    # If we found the final step, check if it's a convergent coupling
    if final_step_node:
        check_convergent_coupling(route)

    return convergent_coupling_found, findings_json
