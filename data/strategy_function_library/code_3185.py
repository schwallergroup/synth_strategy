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


LATE_STAGE_AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with primary amine to imide",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with ammonia to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Acylation of secondary amines with anhydrides",
    "Carboxylic acid to amide conversion",
    "Nitrile to amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving late-stage amide formation. It identifies
    this by checking if a reaction at depth <= 1 produces an amide and matches a name
    from a predefined list of known amide formation reaction types.
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

    late_stage_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_formation, findings_json

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and depth <= 1:  # Check at depths 0 and 1 (final steps)
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Analyzing reaction SMILES: {rsmi}")

                product_smiles = rsmi.split(">")[-1]

                # Check if the product contains an amide group
                has_amide_product = False
                if checker.check_fg("Primary amide", product_smiles):
                    has_amide_product = True
                    if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                if checker.check_fg("Secondary amide", product_smiles):
                    has_amide_product = True
                    if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                if checker.check_fg("Tertiary amide", product_smiles):
                    has_amide_product = True
                    if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                print(f"Product contains amide: {has_amide_product}")

                if has_amide_product:
                    # Check for specific amide formation reactions
                    for reaction_type in LATE_STAGE_AMIDE_FORMATION_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(
                                f"Detected late-stage amide formation as final step: {reaction_type}"
                            )
                            late_stage_amide_formation = True
                            if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            
                            # Add structural constraints if both conditions are met
                            # Positional constraint
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "description": "An amide formation reaction must occur in the last two steps of the synthesis.",
                                    "target_group": "LATE_STAGE_AMIDE_FORMATION_REACTIONS",
                                    "position": {
                                        "variable": "depth",
                                        "operator": "<=",
                                        "value": 1
                                    }
                                }
                            })
                            # Co-occurrence constraint
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "description": "A single reaction step must be one of the specified amide formation reactions AND its product must contain an amide functional group.",
                                    "targets": [
                                        "LATE_STAGE_AMIDE_FORMATION_REACTIONS",
                                        [
                                            "Primary amide",
                                            "Secondary amide",
                                            "Tertiary amide"
                                        ]
                                    ]
                                }
                            })
                            return
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            if (
                not late_stage_amide_formation
            ):  # Stop traversal if we already found what we're looking for
                # New logic for depth calculation
                new_depth = depth
                if node['type'] != 'reaction': # If current node is chemical, depth increases
                    new_depth = depth + 1
                # If current node is reaction, depth remains the same
                dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return late_stage_amide_formation, findings_json
