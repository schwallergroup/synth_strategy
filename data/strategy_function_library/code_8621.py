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


AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Ester with ammonia to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Acyl chloride with ammonia to amide",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route uses a late-stage amide coupling as the final step.
    This looks for a reaction at depth 0 that forms an amide bond.
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
    min_depth_amide_coupling = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal result, min_depth_amide_coupling, findings_json

        print(f"Traversing node type: {node['type']}, depth: {depth}")

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                is_amide_coupling = False
                for reaction_type in AMIDE_COUPLING_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected amide coupling reaction at depth {depth}: {reaction_type}")
                        is_amide_coupling = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        if depth < min_depth_amide_coupling:
                            min_depth_amide_coupling = depth
                        break

                if depth == 0 and is_amide_coupling:
                    result = True

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node['type'] != 'reaction': # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    print("Starting traversal")
    dfs_traverse(route)
    print(f"Min depth amide coupling: {min_depth_amide_coupling}")

    if not result and min_depth_amide_coupling <= 1:
        print(
            f"Found amide coupling at depth {min_depth_amide_coupling}, considering it late-stage"
        )
        result = True
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "amide_coupling",
                "position": "late_stage"
            }
        })

    print(f"Final result: {result}")
    return result, findings_json