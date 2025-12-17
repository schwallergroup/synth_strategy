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
    Checks if the route contains late-stage sulfonamide formation.
    Late-stage means the sulfonamide is formed in the final steps (low depth).
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

    # Track sulfonamide formation reactions and their depths
    sulfonamide_reactions = []
    max_depth = 0
    result = False

    def dfs(node, depth=0):
        nonlocal max_depth, sulfonamide_reactions, findings_json
        max_depth = max(max_depth, depth)

        # Check for sulfonamide formation in reaction nodes
        if node.get("type") == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check if the product contains a sulfonamide group
            if checker.check_fg("Sulfonamide", product):
                if "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")

                # Check if reactants don't have sulfonamide (indicating formation)
                reactants = reactants_part.split(".")
                has_sulfonamide_in_reactants = False

                for reactant in reactants:
                    if checker.check_fg("Sulfonamide", reactant):
                        has_sulfonamide_in_reactants = True
                        break

                # If product has sulfonamide but reactants don't, it's a formation reaction
                if not has_sulfonamide_in_reactants:
                    sulfonamide_reactions.append((depth, rsmi))
                    # Assuming 'Sulfonamide_formation' is the named reaction for this event
                    if "Sulfonamide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Sulfonamide_formation")

        # Recursively check children
        for child in node.get("children", []):
            # New depth calculation logic
            if node.get("type") == "reaction":
                dfs(child, depth)
            else: # Assuming 'chemical' or other types
                dfs(child, depth + 1)

    # Start DFS from the root
    dfs(route)

    # If we found sulfonamide formations, check if any are late-stage
    if sulfonamide_reactions:
        # Sort by depth (ascending)
        sulfonamide_reactions.sort(key=lambda x: x[0])

        # Consider it late-stage if it's in the first third of the synthesis
        late_stage_threshold = max(2, max_depth // 3)

        if sulfonamide_reactions[0][0] <= late_stage_threshold:
            result = True
            # Add the structural constraint if the condition is met
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "Sulfonamide_formation",
                    "position": "late_stage"
                }
            })

    return result, findings_json
