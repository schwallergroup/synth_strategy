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


SUZUKI_REACTION_TYPES = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
    "{Suzuki}",
]

BORYLATION_REACTION_TYPES = [
    "Preparation of boronic acids",
    "Preparation of boronic acids without boronic ether",
    "Preparation of boronic acids from trifluoroborates",
    "Preparation of boronic esters",
    "Synthesis of boronic acids",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy involving a late-stage Suzuki coupling that is preceded by a borylation step. It specifically checks for Suzuki reactions from the list SUZUKI_REACTION_TYPES occurring in the final two steps of the synthesis. It also identifies borylation reactions from the list BORYLATION_REACTION_TYPES occurring at any earlier stage.
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

    # Track if we found the pattern
    found_suzuki = False
    found_borylation = False

    # Track Suzuki coupling and borylation reactions
    suzuki_reactions = []
    borylation_reactions = []

    # Original strategy JSON for structural constraints
    original_strategy_json = {
      "function_id": "code_3258",
      "filepath": "../data/merged_good_perf/code_3258.py",
      "description": "This function detects a synthetic strategy involving a late-stage Suzuki coupling that is preceded by a borylation step. It specifically checks for Suzuki reactions occurring in the final two steps of the synthesis. It also identifies borylation reactions occurring at any earlier stage. The strategy is only confirmed if the Suzuki reaction consumes a reactant containing a boronic acid or boronic ester.",
      "atomic_checks": {
        "named_reactions": [
          "Suzuki coupling with boronic acids",
          "Suzuki coupling with boronic acids OTf",
          "Suzuki coupling with boronic esters",
          "Suzuki coupling with boronic esters OTf",
          "Suzuki coupling with sulfonic esters",
          "{Suzuki}",
          "Preparation of boronic acids",
          "Preparation of boronic acids without boronic ether",
          "Preparation of boronic acids from trifluoroborates",
          "Preparation of boronic esters",
          "Synthesis of boronic acids"
        ],
        "ring_systems": [],
        "functional_groups": [
          "Boronic acid",
          "Boronic ester"
        ]
      },
      "structural_constraints": [
        {
          "type": "co-occurrence",
          "details": {
            "targets": [
              "any_suzuki_reaction",
              "any_borylation_reaction"
            ],
            "description": "The route must contain at least one Suzuki coupling reaction and one borylation reaction."
          }
        },
        {
          "type": "positional",
          "details": {
            "target": "any_suzuki_reaction",
            "position": "last_two_stages",
            "description": "A Suzuki coupling reaction must occur within the last two steps of the synthesis (depth 0 or 1)."
          }
        },
        {
          "type": "sequence",
          "details": {
            "before": "any_borylation_reaction",
            "after": "any_suzuki_reaction",
            "description": "A borylation reaction must occur in an earlier synthetic step (higher depth) than the latest-stage Suzuki coupling."
          }
        },
        {
          "type": "co-occurrence",
          "details": {
            "targets": [
              "any_suzuki_reaction",
              "Boronic acid",
              "Boronic ester"
            ],
            "description": "A Suzuki coupling reaction must have a reactant that contains either a Boronic acid or Boronic ester functional group."
          }
        }
      ]
    }

    def dfs_traverse(node, depth=0):
        nonlocal found_suzuki, found_borylation, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check for Suzuki coupling at depth 0 or 1 (late stage)
                if depth <= 1:
                    is_suzuki = False
                    for rxn_name in SUZUKI_REACTION_TYPES:
                        if checker.check_reaction(rxn_name, rsmi):
                            is_suzuki = True
                            if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                            break

                    if is_suzuki:
                        print(f"Found Suzuki coupling at depth {depth}: {rsmi}")
                        found_suzuki = True
                        suzuki_reactions.append(
                            {"depth": depth, "rsmi": rsmi, "reactants": reactants}
                        )
                        # Record positional constraint if met
                        if depth <= 1:
                            constraint_obj = next((c for c in original_strategy_json["structural_constraints"] if c.get("type") == "positional" and c["details"].get("target") == "any_suzuki_reaction"), None)
                            if constraint_obj and constraint_obj not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(copy.deepcopy(constraint_obj))

                # Check for borylation at any depth
                is_borylation = False
                for rxn_name in BORYLATION_REACTION_TYPES:
                    if checker.check_reaction(rxn_name, rsmi):
                        is_borylation = True
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        break

                if is_borylation:
                    print(f"Found borylation step at depth {depth}: {rsmi}")
                    found_borylation = True
                    borylation_reactions.append(
                        {"depth": depth, "rsmi": rsmi, "product": product_part}
                    )

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)

    # Check if borylation occurs before Suzuki coupling (higher depth)
    borylation_before_suzuki = False

    if found_suzuki and found_borylation:
        # Get minimum depth of Suzuki coupling (latest stage)
        min_suzuki_depth = min(rxn["depth"] for rxn in suzuki_reactions)

        # Check if any borylation occurs at a higher depth (earlier stage)
        for borylation in borylation_reactions:
            if borylation["depth"] > min_suzuki_depth:
                print(
                    f"Confirmed borylation at depth {borylation['depth']} occurs before Suzuki at depth {min_suzuki_depth}"
                )
                borylation_before_suzuki = True
                # Record sequence constraint if met
                constraint_obj = next((c for c in original_strategy_json["structural_constraints"] if c.get("type") == "sequence"), None)
                if constraint_obj and constraint_obj not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(copy.deepcopy(constraint_obj))
                break

    # Check if any Suzuki coupling uses a boronic compound
    suzuki_uses_boronic = False
    if found_suzuki:
        for suzuki in suzuki_reactions:
            found_boronic_acid = False
            found_boronic_ester = False
            for r in suzuki["reactants"]:
                if checker.check_fg("Boronic acid", r):
                    found_boronic_acid = True
                    if "Boronic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Boronic acid")
                if checker.check_fg("Boronic ester", r):
                    found_boronic_ester = True
                    if "Boronic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Boronic ester")
            if found_boronic_acid or found_boronic_ester:
                print(f"Confirmed Suzuki coupling uses a boronic compound")
                suzuki_uses_boronic = True
                # Record co-occurrence constraint if met
                constraint_obj = next((c for c in original_strategy_json["structural_constraints"] if c.get("type") == "co-occurrence" and "Boronic acid" in c["details"].get("targets", [])), None)
                if constraint_obj and constraint_obj not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(copy.deepcopy(constraint_obj))
                break

    # Return True if all conditions are met
    result = found_suzuki and found_borylation and borylation_before_suzuki and suzuki_uses_boronic

    # Record the first co-occurrence constraint if both Suzuki and Borylation are found
    if found_suzuki and found_borylation:
        constraint_obj = next((c for c in original_strategy_json["structural_constraints"] if c.get("type") == "co-occurrence" and "any_suzuki_reaction" in c["details"].get("targets", []) and "any_borylation_reaction" in c["details"].get("targets", [])), None)
        if constraint_obj and constraint_obj not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(copy.deepcopy(constraint_obj))

    print(
        f"Final result: found_suzuki={found_suzuki}, found_borylation={found_borylation}, "
        f"borylation_before_suzuki={borylation_before_suzuki}, suzuki_uses_boronic={suzuki_uses_boronic}, "
        f"returning {result}"
    )
    return result, findings_json