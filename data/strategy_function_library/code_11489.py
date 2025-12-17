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


# Refactoring for Enumeration: Isolate the list of reaction types
AMIDE_FORMATION_REACTION_TYPES = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Acyl chloride with ammonia to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Ester with ammonia to amide",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage amide formation strategy occurring in the final step of the synthesis. 
    The identification is based on the reaction matching one of the specific reaction types 
    defined in the AMIDE_FORMATION_REACTION_TYPES list.
    """
    print(f"Analyzing route for late-stage amide formation strategy")

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Validate route structure
    if not isinstance(route, dict) or "type" not in route:
        print("Invalid route structure")
        return False, findings_json

    has_late_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_amide_formation, findings_json

        print(f"Traversing node at depth {depth}, type: {node.get('type')}")

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate reaction
            print(f"Analyzing late-stage reaction at depth {depth}")

            # Check if metadata and rsmi exist
            if "metadata" not in node or "rsmi" not in node.get("metadata", {}):
                print(f"Missing metadata or rsmi at depth {depth}")
                return

            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Analyzing reaction SMILES: {rsmi}")

                # Check if this is an amide formation reaction using available reaction types
                for reaction_type in AMIDE_FORMATION_REACTION_TYPES:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found late-stage amide formation: {reaction_type}")
                        has_late_amide_formation = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        # Add the structural constraint if this condition is met
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "targets": [
                                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                                    "Carboxylic acid with primary amine to amide",
                                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                                    "Acyl chloride with secondary amine to amide",
                                    "Acyl chloride with ammonia to amide",
                                    "Ester with primary amine to amide",
                                    "Ester with secondary amine to amide",
                                    "Ester with ammonia to amide",
                                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                    "Schotten-Baumann_amide"
                                ],
                                "position": "final_or_penultimate_stage"
                            }
                        })
                        return

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversing the tree
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Late-stage amide formation detected: {has_late_amide_formation}")
    return has_late_amide_formation, findings_json