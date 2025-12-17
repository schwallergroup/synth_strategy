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


LATE_STAGE_COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Heck terminal vinyl",
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Stille reaction_aryl",
    "Negishi coupling",
    "Ullmann condensation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a linear synthesis strategy with a late-stage coupling
    to introduce structural complexity.
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
    
    # Flags to track if specific structural constraints are met
    late_stage_coupling_event_found = False
    sufficient_linear_steps_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_coupling, linear_steps, findings_json, late_stage_coupling_event_found

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")

            # Check if it's a coupling reaction (at late stage)
            if depth <= 2 and len(reactants) >= 2:
                # Check for complex reactants
                complex_reactants = sum(
                    1
                    for r in reactants
                    if Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).GetNumAtoms() > 8
                )

                # Check if it's a known coupling reaction type
                is_coupling_reaction = False
                for rxn_type in LATE_STAGE_COUPLING_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_coupling_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        print(f"Found late-stage coupling reaction: {rxn_type} at depth {depth}")
                        break

                if complex_reactants >= 2 and is_coupling_reaction:
                    late_coupling = True
                    late_stage_coupling_event_found = True
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "any_LATE_STAGE_COUPLING_REACTIONS",
                            "position": "late_stage",
                            "condition": "depth <= 2"
                        }
                    })
                    print(
                        f"Confirmed late-stage coupling at depth {depth} with {complex_reactants} complex reactants"
                    )

            # Count linear steps (single complex reactant)
            if len(reactants) >= 1:
                complex_reactants = sum(
                    1
                    for r in reactants
                    if Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).GetNumAtoms() > 8
                )
                if complex_reactants == 1 and len(reactants) <= 2:
                    linear_steps += 1
                    print(f"Found linear synthesis step at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Summary: late_coupling={late_coupling}, linear_steps={linear_steps}")
    
    result = late_coupling and linear_steps >= 3

    if linear_steps >= 3:
        sufficient_linear_steps_found = True
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "linear_synthesis_step",
                "operator": ">=",
                "value": 3
            }
        })

    if late_stage_coupling_event_found and sufficient_linear_steps_found:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "late_stage_coupling_event",
                    "sufficient_linear_steps"
                ]
            }
        })

    # Return True if we found a late-stage coupling and at least 3 linear steps
    return result, findings_json
