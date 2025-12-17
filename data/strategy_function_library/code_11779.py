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


# Refactoring for Enumeration: Isolate lists of reaction names
ESTERIFICATION_REACTIONS = [
    "Esterification of Carboxylic Acids",
    "Oxidative esterification of primary alcohols",
]

AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Ester with ammonia to amide",
    "Ester with secondary amine to amide",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with secondary amine to amide",
    "Schotten-Baumann_amide",
]

ESTER_DEPROTECTION_REACTIONS = [
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "COOH ethyl deprotection",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route follows a protection-coupling-deprotection strategy,
    specifically looking for esterification followed by amide coupling.
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
    result = False

    # Track reactions in sequence with their associated molecules
    reaction_sequence = []

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for esterification (protection)
                for rxn in ESTERIFICATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        reaction_sequence.append(("esterification", depth, rsmi))
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)

                # Check for amide formation (coupling)
                for rxn in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        reaction_sequence.append(("amide_formation", depth, rsmi))
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)

                # Check for ester hydrolysis (deprotection)
                for rxn in ESTER_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        reaction_sequence.append(("deprotection", depth, rsmi))
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)

            except Exception as e:
                pass

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # Only increase depth when going from chemical to reaction
            next_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Sort reactions by depth to get chronological order (higher depth = earlier stage)
    reaction_sequence.sort(key=lambda x: x[1], reverse=True)
    reaction_types = [r[0] for r in reaction_sequence]

    # Check for protection-coupling-deprotection pattern
    # We need to find all three steps in the correct order
    for i in range(len(reaction_types) - 2):
        if (
            reaction_types[i] == "esterification"
            and reaction_types[i + 1] == "amide_formation"
            and reaction_types[i + 2] == "deprotection"
        ):
            result = True
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "ordered_targets": [
                        "esterification",
                        "amide_formation",
                        "deprotection"
                    ],
                    "description": "Checks for a sequence where esterification is immediately followed by amide formation, which is immediately followed by deprotection in the chronological reaction order."
                }
            })
            # Continue checking for other patterns, but result is already True

    # Also check for protection-coupling pattern (without deprotection)
    for i in range(len(reaction_types) - 1):
        if reaction_types[i] == "esterification" and reaction_types[i + 1] == "amide_formation":
            result = True
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "ordered_targets": [
                        "esterification",
                        "amide_formation"
                    ],
                    "description": "Checks for a sequence where esterification is immediately followed by amide formation in the chronological reaction order."
                }
            })
            # Continue checking for other patterns, but result is already True

    # Check for deprotection-amide_formation pattern (protection might have happened earlier)
    for i in range(len(reaction_types) - 1):
        if reaction_types[i] == "deprotection" and reaction_types[i + 1] == "amide_formation":
            result = True
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "ordered_targets": [
                        "deprotection",
                        "amide_formation"
                    ],
                    "description": "Checks for a sequence where deprotection is immediately followed by amide formation in the chronological reaction order."
                }
            })
            # Continue checking for other patterns, but result is already True

    # Check if we have both esterification and amide_formation in any order
    if "esterification" in reaction_types and "amide_formation" in reaction_types:
        result = True
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "esterification",
                    "amide_formation"
                ],
                "description": "Checks for the presence of both an esterification and an amide formation reaction anywhere in the route, regardless of order."
            }
        })

    return result, findings_json
