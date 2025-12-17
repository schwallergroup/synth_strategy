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
    This function detects a strategy where the final step is an amide formation,
    preceded by an ester hydrolysis, with cyclopropyl fragments incorporated earlier.
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

    # Track if we found the key reactions
    found_amide_formation = False
    found_ester_hydrolysis = False
    found_cyclopropyl_fragments = 0

    # Define structural constraints for later use
    structural_constraints_definitions = [
        {
            "type": "positional",
            "details": {
                "target": "amide_formation",
                "description": "An amide formation reaction must occur within the last three stages of the synthesis (depth <= 2).",
                "max_depth": 2
            }
        },
        {
            "type": "positional",
            "details": {
                "target": "ester_hydrolysis",
                "description": "An ester hydrolysis reaction must occur at an intermediate stage of the synthesis (depth between 1 and 3).",
                "min_depth": 1,
                "max_depth": 3
            }
        },
        {
            "type": "count",
            "details": {
                "target": "cyclopropane",
                "operator": ">",
                "value": 0,
                "description": "At least one cyclopropane ring must be present in a molecule or reaction within the route."
            }
        },
        {
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "amide_formation",
                    "ester_hydrolysis",
                    "cyclopropane"
                ],
                "description": "The overall strategy requires the successful identification of a late-stage amide formation, an intermediate-stage ester hydrolysis, and the presence of at least one cyclopropane fragment."
            }
        }
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_formation, found_ester_hydrolysis, found_cyclopropyl_fragments, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide formation at late stage (depth <= 2)
                if depth <= 2:
                    # Check for amide formation reaction directly
                    amide_formation_reactions = [
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                        "Carboxylic acid with primary amine to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Acyl chloride with secondary amine to amide",
                        "Ester with primary amine to amide",
                        "Ester with secondary amine to amide",
                        "Acyl chloride with ammonia to amide",
                        "Ester with ammonia to amide",
                    ]

                    for reaction_type in amide_formation_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            found_amide_formation = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            print(f"Found amide formation at depth {depth}: {reaction_type}")
                            break

                # Check for ester hydrolysis at intermediate stage (1 <= depth <= 3)
                if 1 <= depth <= 3:
                    # Check for ester hydrolysis reaction directly
                    hydrolysis_reactions = [
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                        "Ester saponification (methyl deprotection)",
                        "Ester saponification (alkyl deprotection)",
                        "COOH ethyl deprotection",
                        "Deprotection of carboxylic acid",
                    ]

                    for reaction_type in hydrolysis_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            found_ester_hydrolysis = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            print(f"Found ester hydrolysis at depth {depth}: {reaction_type}")
                            break

                # Check for cyclopropyl fragments at any step
                cyclopropyl_found_in_reaction = False
                for reactant in reactants:
                    if checker.check_ring("cyclopropane", reactant):
                        found_cyclopropyl_fragments += 1
                        if "cyclopropane" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("cyclopropane")
                        cyclopropyl_found_in_reaction = True
                if checker.check_ring("cyclopropane", product):
                    found_cyclopropyl_fragments += 1
                    if "cyclopropane" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("cyclopropane")
                    cyclopropyl_found_in_reaction = True

                if cyclopropyl_found_in_reaction:
                    print(f"Found cyclopropyl fragment at depth {depth}")

        # Check for cyclopropyl fragments in molecule nodes too
        elif node["type"] == "mol" and node["smiles"]:
            if checker.check_ring("cyclopropane", node["smiles"]):
                found_cyclopropyl_fragments += 1
                if "cyclopropane" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("cyclopropane")
                print(f"Found cyclopropyl fragment in molecule at depth {depth}")

        # Recursively process children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not 'reaction' (e.g., 'chemical')
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # The strategy is present if we found amide formation at the end,
    # ester hydrolysis before that, and at least one cyclopropyl fragment
    strategy_present = (
        found_amide_formation and found_ester_hydrolysis and found_cyclopropyl_fragments > 0
    )

    print(f"Amide formation: {found_amide_formation}")
    print(f"Ester hydrolysis: {found_ester_hydrolysis}")
    print(f"Cyclopropyl fragments: {found_cyclopropyl_fragments}")

    # Record structural constraints if met
    if found_amide_formation:
        for constraint in structural_constraints_definitions:
            if constraint["type"] == "positional" and constraint["details"]["target"] == "amide_formation":
                findings_json["structural_constraints"].append(constraint)
                break
    if found_ester_hydrolysis:
        for constraint in structural_constraints_definitions:
            if constraint["type"] == "positional" and constraint["details"]["target"] == "ester_hydrolysis":
                findings_json["structural_constraints"].append(constraint)
                break
    if found_cyclopropyl_fragments > 0:
        for constraint in structural_constraints_definitions:
            if constraint["type"] == "count" and constraint["details"]["target"] == "cyclopropane":
                findings_json["structural_constraints"].append(constraint)
                break
    if strategy_present:
        for constraint in structural_constraints_definitions:
            if constraint["type"] == "co-occurrence" and "amide_formation" in constraint["details"]["targets"] and "ester_hydrolysis" in constraint["details"]["targets"] and "cyclopropane" in constraint["details"]["targets"]:
                findings_json["structural_constraints"].append(constraint)
                break

    return strategy_present, findings_json
