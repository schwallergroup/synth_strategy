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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

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


NAMED_CROSS_COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters OTf",
    "Stille reaction_aryl",
    "Stille reaction_vinyl",
    "Negishi coupling",
    "Heck terminal vinyl",
    "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl halide",
    "Buchwald-Hartwig",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Ullmann condensation",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Aryllithium cross-coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a known, named cross-coupling reaction occurs in the late stages of a synthesis
    (defined as the final third of the steps). This strategy is indicative of a convergent
    synthesis plan where major fragments are joined near the end. The specific reactions
    checked are defined in the NAMED_CROSS_COUPLING_REACTIONS list.
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

    fragment_coupling_detected = False
    total_depth = 0

    # First pass to determine total depth
    def get_max_depth(node, current_depth=0):
        nonlocal total_depth
        total_depth = max(total_depth, current_depth)
        for child in node.get("children", []):
            get_max_depth(child, current_depth + 1)

    get_max_depth(route)

    # Define what counts as "late stage" based on total synthesis length
    late_stage_threshold = total_depth // 3 if total_depth > 3 else 1

    def dfs_traverse(node, depth=0):
        nonlocal fragment_coupling_detected, findings_json

        if node["type"] == "reaction" and depth <= late_stage_threshold:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                # Skip if we don't have at least 2 reactants
                if len(reactants) < 2:
                    for child in node.get("children", []):
                        # Depth remains the same when going from reaction to chemical
                        dfs_traverse(child, depth)
                    return

                # Check if this is a known coupling reaction
                is_coupling_reaction_found = False
                for name in NAMED_CROSS_COUPLING_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        is_coupling_reaction_found = True
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                        break

                if is_coupling_reaction_found:
                    fragment_coupling_detected = True
                    # Add the structural constraint if a late-stage coupling is detected
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "any_named_cross_coupling",
                            "position": "late_stage"
                        }
                    })
                    return

        for child in node.get("children", []):
            if not fragment_coupling_detected:  # Stop traversal if we already found a coupling
                # New logic: depth increases only from chemical to reaction
                # If current node is 'reaction', depth remains the same for children (chemicals)
                # If current node is 'chemical', depth increases for children (reactions)
                next_depth = depth + 1 if node["type"] != "reaction" else depth
                dfs_traverse(child, next_depth)

    dfs_traverse(route)
    return fragment_coupling_detected, findings_json
