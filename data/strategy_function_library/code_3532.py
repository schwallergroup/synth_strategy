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
    This function detects a pattern of esterification-hydrolysis-esterification
    sequences in the synthetic route.
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

    # Track reactions by depth
    reactions_by_depth = {}

    def dfs_traverse(node, current_depth=0):
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            try:
                depth = int(
                    re.search(
                        r"Depth: (\d+)", node["metadata"].get("ID", f"Depth: {current_depth}")
                    ).group(1)
                )
            except (AttributeError, ValueError):
                depth = current_depth

            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactions_by_depth[depth] = rsmi
            print(f"Found reaction at depth {depth}: {rsmi}")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for its children (which are chemicals)
                dfs_traverse(child, current_depth)
            else:
                # If current node is a chemical, depth increases for its children (which are reactions)
                dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Sort reactions by depth (reverse order for retrosynthetic analysis)
    sorted_depths = sorted(reactions_by_depth.keys(), reverse=True)
    print(f"Sorted depths (reverse order for retrosynthesis): {sorted_depths}")

    # Check the sequence of reactions
    reaction_sequence = []
    for depth in sorted_depths:
        rsmi = reactions_by_depth[depth]
        print(f"Checking reaction at depth {depth}: {rsmi}")

        # Split reaction SMILES to get reactants and products
        try:
            # Check for esterification (carboxylic acid → ester)
            if checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                reaction_sequence.append("esterification")
                if "Esterification of Carboxylic Acids" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Esterification of Carboxylic Acids")
                print(f"Detected esterification at depth {depth}")

            # Check for hydrolysis (ester → carboxylic acid)
            elif checker.check_reaction(
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
            ):
                reaction_sequence.append("hydrolysis")
                if "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters")
                print(f"Detected hydrolysis at depth {depth}")

        except Exception as e:
            print(f"Error processing reaction at depth {depth}: {e}")

    print(f"Final reaction sequence: {reaction_sequence}")

    # Check if the sequence contains esterification-hydrolysis-esterification
    pattern_detected = False

    # Check for esterification-hydrolysis-esterification pattern
    for i in range(len(reaction_sequence) - 2):
        if (
            reaction_sequence[i] == "esterification"
            and reaction_sequence[i + 1] == "hydrolysis"
            and reaction_sequence[i + 2] == "esterification"
        ):
            pattern_detected = True
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "ordered_events": [
                        "Esterification of Carboxylic Acids",
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                        "Esterification of Carboxylic Acids"
                    ],
                    "description": "Checks for a sequence of esterification, followed by hydrolysis, followed by another esterification in the retrosynthetic path."
                }
            })
            print(
                f"Esterification-hydrolysis-esterification pattern detected at positions {i}, {i+1}, {i+2}"
            )
            break

    return pattern_detected, findings_json