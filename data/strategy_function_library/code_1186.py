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


PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
    "Alcohol protection with silyl ethers",
    "Protection of carboxylic acid",
]

COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Heck terminal vinyl",
    "Negishi coupling",
    "Stille reaction_aryl",
    "Ullmann condensation",
    "Goldberg coupling",
    "{Suzuki}",
    "{N-arylation_heterocycles}",
    "{Buchwald-Hartwig}",
]

DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Alcohol deprotection from silyl ethers",
    "Alcohol deprotection from silyl ethers (double)",
    "Alcohol deprotection from silyl ethers (diol)",
    "Deprotection of carboxylic acid",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "Hydroxyl benzyl deprotection",
    "Carboxyl benzyl deprotection",
    "COOH ethyl deprotection",
    "Phthalimide deprotection",
    "Tert-butyl deprotection of amine",
    "TMS deprotection from alkyne",
    "N-glutarimide deprotection",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthesis route follows a protection-coupling-deprotection pattern. This is determined by identifying reactions from predefined lists of named protection, coupling, and deprotection reactions and checking if they occur in the correct synthetic order (protection -> coupling -> deprotection).
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

    # Track protection, coupling, and deprotection events in order
    events = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for protection event
                for protection_rxn in PROTECTION_REACTIONS:
                    if checker.check_reaction(protection_rxn, rsmi):
                        events.append(("protection", depth, rsmi))
                        findings_json["atomic_checks"]["named_reactions"].append(protection_rxn)
                        print(f"Found protection at depth {depth}: {protection_rxn}")
                        break

                # Check for coupling event
                for coupling_rxn in COUPLING_REACTIONS:
                    if checker.check_reaction(coupling_rxn, rsmi):
                        events.append(("coupling", depth, rsmi))
                        findings_json["atomic_checks"]["named_reactions"].append(coupling_rxn)
                        print(f"Found coupling at depth {depth}: {coupling_rxn}")
                        break

                # Check for deprotection event
                for deprotection_rxn in DEPROTEPTION_REACTIONS:
                    if checker.check_reaction(deprotection_rxn, rsmi):
                        events.append(("deprotection", depth, rsmi))
                        findings_json["atomic_checks"]["named_reactions"].append(deprotection_rxn)
                        print(f"Found deprotection at depth {depth}: {deprotection_rxn}")
                        break
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a chemical node
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    dfs_traverse(route)

    # Check for protection-coupling-deprotection pattern
    # Sort events by depth to ensure chronological order in the synthesis direction
    # Note: In retrosynthesis, higher depth = earlier in actual synthesis
    events.sort(key=lambda x: x[1], reverse=True)

    print(f"Sorted events: {events}")

    result = False

    if not events:
        print("No protection, coupling, or deprotection events found")
        return result, findings_json

    # Extract depths for each event type
    protection_depths = [d for t, d, _ in events if t == "protection"]
    coupling_depths = [d for t, d, _ in events if t == "coupling"]
    deprotection_depths = [d for t, d, _ in events if t == "deprotection"]

    # Check if we have all three event types
    if protection_depths and coupling_depths and deprotection_depths:
        # Add co-occurrence constraint if all types are found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "protection_reaction",
                    "coupling_reaction",
                    "deprotection_reaction"
                ],
                "description": "The route must contain at least one reaction from each of the protection, coupling, and deprotection categories."
            }
        })

        # Check if they occur in the correct order in the synthesis direction
        # In synthesis direction (reverse of retrosynthesis), we want:
        # protection (highest depth) -> coupling -> deprotection (lowest depth)
        if max(protection_depths) > max(coupling_depths) > max(deprotection_depths):
            print("Found protection-coupling-deprotection pattern in correct order")
            result = True
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "ordered_events": [
                        "protection_reaction",
                        "coupling_reaction",
                        "deprotection_reaction"
                    ],
                    "description": "A protection reaction must occur synthetically before a coupling reaction, which must occur before a deprotection reaction. The check compares the latest (in synthesis time) occurrence of each reaction type."
                }
            })

    # Alternative check: look for the pattern in the sequence of events
    found_protection = False
    found_coupling = False
    found_deprotection = False

    for event_type, depth, _ in events:
        if event_type == "protection" and not found_protection:
            found_protection = True
        elif event_type == "coupling" and found_protection and not found_coupling:
            found_coupling = True
        elif event_type == "deprotection" and found_coupling and not found_deprotection:
            found_deprotection = True

    if found_protection and found_coupling and found_deprotection:
        print("Found protection-coupling-deprotection pattern in sequential order")
        # This condition also implies the co-occurrence constraint is met
        if {
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "protection_reaction",
                    "coupling_reaction",
                    "deprotection_reaction"
                ],
                "description": "The route must contain at least one reaction from each of the protection, coupling, and deprotection categories."
            }
        } not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "protection_reaction",
                        "coupling_reaction",
                        "deprotection_reaction"
                    ],
                    "description": "The route must contain at least one reaction from each of the protection, coupling, and deprotection categories."
                }
            })
        
        # This condition also implies the sequence constraint is met
        if {
            "type": "sequence",
            "details": {
                "ordered_events": [
                    "protection_reaction",
                    "coupling_reaction",
                    "deprotection_reaction"
                ],
                "description": "A protection reaction must occur synthetically before a coupling reaction, which must occur before a deprotection reaction. The check compares the latest (in synthesis time) occurrence of each reaction type."
            }
        } not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "ordered_events": [
                        "protection_reaction",
                        "coupling_reaction",
                        "deprotection_reaction"
                    ],
                    "description": "A protection reaction must occur synthetically before a coupling reaction, which must occur before a deprotection reaction. The check compares the latest (in synthesis time) occurrence of each reaction type."
                }
            })
        result = True

    return result, findings_json