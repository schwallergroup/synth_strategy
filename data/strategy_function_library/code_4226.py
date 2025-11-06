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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a sequential nitrogen functionalization strategy:
    nitro → primary amine → secondary amine → tertiary amide/urea
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
    found_nitro_reduction = False
    found_reductive_amination = False
    found_amide_formation = False

    # Track the depth at which each transformation occurs
    nitro_reduction_depth = -1
    reductive_amination_depth = -1
    amide_formation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_nitro_reduction, found_reductive_amination, found_amide_formation
        nonlocal nitro_reduction_depth, reductive_amination_depth, amide_formation_depth
        nonlocal findings_json

        if node["type"] == "reaction":
            try:
                # Extract rsmi
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for nitro reduction
                nitro_reduction_reactions = ["Reduction of nitro groups to amines"]
                for r_name in nitro_reduction_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        found_nitro_reduction = True
                        nitro_reduction_depth = depth
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                # Check for reductive amination
                reductive_amination_reactions = [
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "reductive amination"
                ]
                for r_name in reductive_amination_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        found_reductive_amination = True
                        reductive_amination_depth = depth
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                # Check for amide/urea formation from a secondary amine
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with secondary amine to amide",
                    "Acylation of secondary amines",
                    "Acylation of secondary amines with anhydrides",
                    "Schotten-Baumann_amide",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids"
                ]
                urea_formation_reactions = [
                    "Urea synthesis via isocyanate and secondary amine",
                    "urea"
                ]

                is_amide_formation = False
                for r_name in amide_formation_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        is_amide_formation = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                is_urea_formation = False
                for r_name in urea_formation_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        is_urea_formation = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                if is_amide_formation or is_urea_formation:
                    found_amide_formation = True
                    amide_formation_depth = depth

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a chemical node
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    result = False
    # Check if we found the complete pattern in the correct order
    if found_nitro_reduction and found_reductive_amination and found_amide_formation:
        # Check if the transformations are in the correct order (higher depth = earlier in synthesis)
        if nitro_reduction_depth > reductive_amination_depth > amide_formation_depth:
            result = True
            # Add structural constraints if the full pattern is found
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "nitro_reduction",
                        "reductive_amination",
                        "amide_or_urea_formation"
                    ],
                    "description": "The route must contain at least one reaction from each of the three conceptual categories: nitro reduction, reductive amination, and amide/urea formation. The specific reactions for each category are listed in atomic_checks."
                }
            })
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "ordered_events": [
                        "nitro_reduction",
                        "reductive_amination",
                        "amide_or_urea_formation"
                    ],
                    "description": "The identified reactions must occur in a specific chronological order during the synthesis: nitro reduction must precede reductive amination, which must precede amide or urea formation."
                }
            })

    return result, findings_json
