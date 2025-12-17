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


OXIDATION_REACTIONS = [
    "Oxidation of aldehydes to carboxylic acids",
    "Oxidation of ketone to carboxylic acid",
    "Oxidation of alcohol to carboxylic acid",
    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
    "Oxidation of alkene to carboxylic acid",
    "Oxidation of alkene to aldehyde",
    "Oxidative esterification of primary alcohols",
    "Oxidation of alcohol and aldehyde to ester",
    "Quinone formation",
    "Oxidation of boronic acids",
    "Oxidation of boronic esters",
]

REDUCTION_REACTIONS = [
    "Reduction of aldehydes and ketones to alcohols",
    "Reduction of ester to primary alcohol",
    "Reduction of ketone to secondary alcohol",
    "Reduction of carboxylic acid to primary alcohol",
    "Reduction of nitrile to amine",
    "Reduction of primary amides to amines",
    "Reduction of secondary amides to amines",
    "Reduction of tertiary amides to amines",
    "Reduction of nitro groups to amines",
    "Hydrogenation (double to single)",
    "Hydrogenation (triple to double)",
    "Arene hydrogenation",
    "Azide to amine reduction (Staudinger)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis uses an early oxidation followed by later reductions.
    In retrosynthetic planning, early stage = high depth, late stage = low depth.
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

    oxidation_depth = -1
    reduction_depths = []
    result = False

    def dfs_traverse(node, current_depth=0):
        nonlocal oxidation_depth, reduction_depths, findings_json

        # Try to get depth from metadata
        depth = current_depth
        if node.get("metadata", {}).get("depth", None) is not None:
            depth = node["metadata"]["depth"]
        elif "metadata" in node and "ID" in node["metadata"]:
            depth_info = str(node["metadata"]["ID"])
            if "Depth:" in depth_info:
                try:
                    depth = int(depth_info.split("Depth:")[1].split()[0])
                except Exception as e:
                    print(f"Error extracting depth from ID: {e}")

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            try:
                # Check for oxidation reactions
                is_oxidation = False
                for rxn_type in OXIDATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_oxidation = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        print(f"Detected oxidation reaction '{rxn_type}' at depth {depth}")
                        break

                if is_oxidation and (oxidation_depth == -1 or depth > oxidation_depth):
                    oxidation_depth = depth
                    print(f"Updated earliest oxidation depth to {depth}")

                # Check for reduction reactions
                is_reduction = False
                for rxn_type in REDUCTION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_reduction = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        print(f"Detected reduction reaction '{rxn_type}' at depth {depth}")
                        break

                if is_reduction:
                    reduction_depths.append(depth)
                    print(f"Added reduction at depth {depth}")

            except Exception as e:
                print(f"Error processing SMILES in early_stage_redox_sequence: {e}")

        # Continue traversing
        for child in node.get("children", []):
            # New logic: depth only increases when traversing from chemical to reaction
            # Depth remains the same when traversing from reaction to chemical
            next_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if oxidation occurs early and reduction occurs later
    if oxidation_depth != -1 and reduction_depths:
        # In retrosynthetic direction, higher depth = earlier in synthesis
        # So oxidation should have higher depth (earlier) than reduction (later)
        min_reduction_depth = min(reduction_depths)
        if oxidation_depth > min_reduction_depth:
            print(
                f"Detected early oxidation (depth {oxidation_depth}) followed by later reduction (depth {min_reduction_depth})"
            )
            result = True
            # Add structural constraints if the condition is met
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "description": "The route must contain at least one oxidation reaction and at least one reduction reaction.",
                    "target_groups": [
                        {
                            "name": "oxidation_reactions",
                            "members": [
                                "Oxidation of aldehydes to carboxylic acids",
                                "Oxidation of ketone to carboxylic acid",
                                "Oxidation of alcohol to carboxylic acid",
                                "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                                "Oxidation of alkene to carboxylic acid",
                                "Oxidation of alkene to aldehyde",
                                "Oxidative esterification of primary alcohols",
                                "Oxidation of alcohol and aldehyde to ester",
                                "Quinone formation",
                                "Oxidation of boronic acids",
                                "Oxidation of boronic esters"
                            ],
                            "min_count": 1
                        },
                        {
                            "name": "reduction_reactions",
                            "members": [
                                "Reduction of aldehydes and ketones to alcohols",
                                "Reduction of ester to primary alcohol",
                                "Reduction of ketone to secondary alcohol",
                                "Reduction of carboxylic acid to primary alcohol",
                                "Reduction of nitrile to amine",
                                "Reduction of primary amides to amines",
                                "Reduction of secondary amides to amines",
                                "Reduction of tertiary amides to amines",
                                "Reduction of nitro groups to amines",
                                "Hydrogenation (double to single)",
                                "Hydrogenation (triple to double)",
                                "Arene hydrogenation",
                                "Azide to amine reduction (Staudinger)"
                            ],
                            "min_count": 1
                        }
                    ]
                }
            })
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": {
                        "target_group": "oxidation_reactions"
                    },
                    "after": {
                        "target_group": "reduction_reactions"
                    },
                    "condition": "The maximum depth of any reaction in the 'before' group must be greater than the minimum depth of any reaction in the 'after' group."
                }
            })
        else:
            print(
                f"Found oxidation and reduction, but oxidation (depth {oxidation_depth}) is not earlier than reduction (depth {min_reduction_depth})"
            )
    else:
        if oxidation_depth == -1:
            print("No oxidation reactions detected")
        if not reduction_depths:
            print("No reduction reactions detected")

    return result, findings_json
