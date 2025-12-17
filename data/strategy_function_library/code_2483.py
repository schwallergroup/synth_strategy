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
    This function detects a strategy of sequential functionalization of an indole scaffold
    with preservation of the core structure throughout the synthesis.
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

    # Track key features
    indole_present = False
    sequential_functionalization = 0
    functionalization_types = set()

    def dfs_traverse(node, depth=0):
        nonlocal indole_present, sequential_functionalization, functionalization_types, findings_json

        if node["type"] == "mol":
            # Check if molecule contains indole core
            if checker.check_ring("indole", node["smiles"]):
                indole_present = True
                findings_json["atomic_checks"]["ring_systems"].append("indole")
                print(f"Depth {depth}: Indole detected in molecule: {node['smiles']}")

        elif node["type"] == "reaction" and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if indole is preserved in the reaction
            product_has_indole = checker.check_ring("indole", product)
            reactants_have_indole = any(checker.check_ring("indole", r) for r in reactants)

            if product_has_indole and reactants_have_indole:
                print(f"Depth {depth}: Indole preserved in reaction: {rsmi}")

                # Check for various functionalization reactions
                # Nitration
                nitration_reactions = [
                    "Aromatic nitration with HNO3",
                    "Aromatic nitration with NO3 salt",
                    "Aromatic nitration with NO2 salt",
                    "Aromatic nitration with alkyl NO2"
                ]
                nitration_detected = False
                for r_name in nitration_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        sequential_functionalization += 1
                        functionalization_types.add("nitration")
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        nitration_detected = True
                        print(f"Depth {depth}: Detected nitration reaction")
                        break

                # Halogenation
                halogenation_reactions = [
                    "Aromatic fluorination",
                    "Aromatic chlorination",
                    "Aromatic bromination",
                    "Aromatic iodination"
                ]
                halogenation_detected = False
                if not nitration_detected: # Only check if nitration wasn't detected
                    for r_name in halogenation_reactions:
                        if checker.check_reaction(r_name, rsmi):
                            sequential_functionalization += 1
                            functionalization_types.add("halogenation")
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                            halogenation_detected = True
                            print(f"Depth {depth}: Detected halogenation reaction")
                            break

                # Cyanation (check for nitrile introduction)
                if not nitration_detected and not halogenation_detected: # Only check if previous weren't detected
                    if checker.check_fg("Nitrile", product) and not any(
                        checker.check_fg("Nitrile", r) for r in reactants
                    ):
                        sequential_functionalization += 1
                        functionalization_types.add("cyanation")
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                        findings_json["atomic_checks"]["named_reactions"].append("functional_group_formation") # Generic reaction for FG formation
                        print(f"Depth {depth}: Detected cyano group introduction")

                # Acylation
                if not nitration_detected and not halogenation_detected and not (checker.check_fg("Nitrile", product) and not any(checker.check_fg("Nitrile", r) for r in reactants)):
                    if checker.check_reaction("Friedel-Crafts acylation", rsmi):
                        sequential_functionalization += 1
                        functionalization_types.add("acylation")
                        findings_json["atomic_checks"]["named_reactions"].append("Friedel-Crafts acylation")
                        print(f"Depth {depth}: Detected acylation reaction")

                # Alkylation
                if not nitration_detected and not halogenation_detected and not (checker.check_fg("Nitrile", product) and not any(checker.check_fg("Nitrile", r) for r in reactants)) and not checker.check_reaction("Friedel-Crafts acylation", rsmi):
                    if checker.check_reaction("Friedel-Crafts alkylation", rsmi):
                        sequential_functionalization += 1
                        functionalization_types.add("alkylation")
                        findings_json["atomic_checks"]["named_reactions"].append("Friedel-Crafts alkylation")
                        print(f"Depth {depth}: Detected alkylation reaction")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., 'mol'), increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Determine if strategy is present
    strategy_present = (
        indole_present
        and sequential_functionalization >= 2
        and len(functionalization_types) >= 2
    )

    print(f"Indole present: {indole_present}")
    print(f"Sequential functionalizations: {sequential_functionalization}")
    print(f"Functionalization types: {functionalization_types}")
    print(f"Strategy detected: {strategy_present}")

    # Populate structural constraints based on final flags
    if indole_present:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "indole"
                ],
                "comment": "The route must contain an indole scaffold."
            }
        })
    if sequential_functionalization >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "indole_preserving_functionalization",
                "operator": ">=",
                "value": 2,
                "comment": "Counts the number of functionalization reactions (nitration, halogenation, cyanation, acylation, alkylation) that preserve the indole core."
            }
        })
    if len(functionalization_types) >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "unique_functionalization_type",
                "operator": ">=",
                "value": 2,
                "comment": "Counts the number of distinct categories of functionalization reactions."
            }
        })

    return strategy_present, findings_json
