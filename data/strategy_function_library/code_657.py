from typing import Tuple, Dict, List
import copy
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
    This function detects if the synthesis involves building a quinazoline scaffold
    followed by functional group modifications and late-stage coupling.
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

    # Track the key features of this strategy
    has_quinazoline_formation = False
    has_functional_group_modification = False
    has_late_stage_coupling = False

    # Original strategy JSON for structural constraints lookup
    original_strategy_json = {
      "function_id": "code_657",
      "filepath": "../data/merged_good_perf/code_657.py",
      "description": "This function detects if the synthesis involves building a quinazoline scaffold followed by functional group modifications and late-stage coupling. The strategy is considered successful if at least two of these three features are detected.",
      "atomic_checks": {
        "named_reactions": [
          "ring_formation",
          "{Niementowski_quinazoline}",
          "Alcohol to chloride_SOCl2",
          "Alcohol to chloride_HCl",
          "Alcohol to chloride_Other",
          "Appel reaction",
          "Nitrile to amide",
          "Oxidation of aldehydes to carboxylic acids",
          "Reduction of nitro groups to amines",
          "Esterification of Carboxylic Acids",
          "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
          "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
          "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
          "{Buchwald-Hartwig}",
          "Goldberg coupling"
        ],
        "ring_systems": [
          "quinazoline"
        ],
        "functional_groups": []
      },
      "structural_constraints": [
        {
          "type": "positional",
          "details": {
            "description": "A quinazoline ring must be formed (present in product, absent in reactants) in an early stage of the synthesis.",
            "target": "ring_formation:quinazoline",
            "position": "depth >= 2"
          }
        },
        {
          "type": "positional",
          "details": {
            "description": "A functional group modification must occur in a middle stage of the synthesis.",
            "target_group": [
              "Alcohol to chloride_SOCl2",
              "Alcohol to chloride_HCl",
              "Alcohol to chloride_Other",
              "Appel reaction",
              "Nitrile to amide",
              "Oxidation of aldehydes to carboxylic acids",
              "Reduction of nitro groups to amines",
              "Esterification of Carboxylic Acids"
            ],
            "position": "1 <= depth < 3"
          }
        },
        {
          "type": "positional",
          "details": {
            "description": "A C-N coupling reaction must occur in a late stage of the synthesis.",
            "target_group": [
              "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
              "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
              "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
              "{Buchwald-Hartwig}",
              "Goldberg coupling"
            ],
            "position": "depth <= 1"
          }
        },
        {
          "type": "count",
          "details": {
            "description": "The overall strategy is valid if at least two of the three key features (quinazoline formation, functional group modification, late-stage coupling) are present.",
            "target": [
              "has_quinazoline_formation",
              "has_functional_group_modification",
              "has_late_stage_coupling"
            ],
            "operator": ">=",
            "value": 2
          }
        }
      ]
    }

    def dfs_traverse(node, depth=0):
        nonlocal has_quinazoline_formation, has_functional_group_modification, has_late_stage_coupling, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for quinazoline formation (early stage)
                if depth >= 2:  # Early stage
                    # Check if product has quinazoline
                    if checker.check_ring("quinazoline", product_smiles):
                        findings_json["atomic_checks"]["ring_systems"].append("quinazoline")
                        # Check if reactants don't have quinazoline
                        reactants_have_quinazoline = False
                        for reactant_smiles in reactants_smiles:
                            if checker.check_ring("quinazoline", reactant_smiles):
                                reactants_have_quinazoline = True
                                break

                        if not reactants_have_quinazoline:
                            print(f"Quinazoline formation detected at depth {depth}")
                            # Check if it's specifically a Niementowski quinazoline formation
                            if checker.check_reaction("{Niementowski_quinazoline}", rsmi):
                                findings_json["atomic_checks"]["named_reactions"].append("{Niementowski_quinazoline}")
                                print(
                                    f"Niementowski quinazoline formation confirmed at depth {depth}"
                                )
                                has_quinazoline_formation = True
                            else:
                                # Alternative check for quinazoline formation
                                has_quinazoline_formation = True
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation") # Generic ring formation
                                print(f"General quinazoline formation detected at depth {depth}")
                            
                            # Add structural constraint for quinazoline formation
                            if has_quinazoline_formation:
                                for constraint in original_strategy_json["structural_constraints"]:
                                    if constraint["type"] == "positional" and constraint["details"]["target"] == "ring_formation:quinazoline":
                                        findings_json["structural_constraints"].append(copy.deepcopy(constraint))
                                        break

                # Check for functional group modifications (middle stage)
                if 1 <= depth < 3:  # Middle stage
                    fg_mod_detected = False
                    # Check for alcohol to halide conversion
                    alcohol_to_halide_reactions = [
                        "Alcohol to chloride_SOCl2",
                        "Alcohol to chloride_HCl",
                        "Alcohol to chloride_Other",
                        "Appel reaction"
                    ]
                    for r_name in alcohol_to_halide_reactions:
                        if checker.check_reaction(r_name, rsmi):
                            print(f"Alcohol to halide conversion detected at depth {depth}")
                            has_functional_group_modification = True
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                            fg_mod_detected = True
                            break

                    # Check for nitrile to amide conversion
                    if checker.check_reaction("Nitrile to amide", rsmi):
                        print(f"Nitrile to amide conversion detected at depth {depth}")
                        has_functional_group_modification = True
                        findings_json["atomic_checks"]["named_reactions"].append("Nitrile to amide")
                        fg_mod_detected = True

                    # Check for other common functional group modifications
                    other_fg_mod_reactions = [
                        "Oxidation of aldehydes to carboxylic acids",
                        "Reduction of nitro groups to amines",
                        "Esterification of Carboxylic Acids"
                    ]
                    for r_name in other_fg_mod_reactions:
                        if checker.check_reaction(r_name, rsmi):
                            print(f"Other functional group modification detected at depth {depth}")
                            has_functional_group_modification = True
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                            fg_mod_detected = True
                            break
                    
                    if fg_mod_detected:
                        for constraint in original_strategy_json["structural_constraints"]:
                            if constraint["type"] == "positional" and "target_group" in constraint["details"] and constraint["details"]["position"] == "1 <= depth < 3":
                                findings_json["structural_constraints"].append(copy.deepcopy(constraint))
                                break

                # Check for late-stage coupling
                if depth <= 1:  # Late stage
                    late_stage_coupling_reactions = [
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        "{Buchwald-Hartwig}",
                        "Goldberg coupling"
                    ]
                    for r_name in late_stage_coupling_reactions:
                        if checker.check_reaction(r_name, rsmi):
                            print(f"Late-stage aromatic amine coupling detected at depth {depth}")
                            has_late_stage_coupling = True
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                            
                            for constraint in original_strategy_json["structural_constraints"]:
                                if constraint["type"] == "positional" and "target_group" in constraint["details"] and constraint["details"]["position"] == "depth <= 1":
                                    findings_json["structural_constraints"].append(copy.deepcopy(constraint))
                                    break
                            break

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is chemical, depth increases
            next_depth = depth + 1

        # Traverse children with the determined next_depth
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # The strategy is present if all three key features are detected
    print(f"Quinazoline formation: {has_quinazoline_formation}")
    print(f"Functional group modification: {has_functional_group_modification}")
    print(f"Late-stage coupling: {has_late_stage_coupling}")

    # Check if any two of the three features are detected (relaxed condition)
    features_count = sum(
        [has_quinazoline_formation, has_functional_group_modification, has_late_stage_coupling]
    )
    result = features_count >= 2  # Return true if at least 2 features are detected

    if result:
        # Add the count structural constraint if the overall condition is met
        for constraint in original_strategy_json["structural_constraints"]:
            if constraint["type"] == "count":
                findings_json["structural_constraints"].append(copy.deepcopy(constraint))
                break

    return result, findings_json
