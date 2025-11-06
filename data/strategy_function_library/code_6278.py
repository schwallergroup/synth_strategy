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
    Detects if the final step of the synthesis is a sulfonamide formation
    via a Schotten-Baumann type reaction between a sulfonyl chloride and an amine.
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

    final_step_is_sulfonamide = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_sulfonamide, findings_json

        if node["type"] == "reaction" and depth == 0:
            try:
                # Extract reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is a sulfonamide formation reaction
                if checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                ):
                    print(f"Detected sulfonamide formation at depth {depth}")
                    final_step_is_sulfonamide = True
                    findings_json["atomic_checks"]["named_reactions"].append("Sulfonamide synthesis (Schotten-Baumann) primary amine")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                            "position": "last_stage"
                        }
                    })
                elif checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                ):
                    print(f"Detected sulfonamide formation at depth {depth}")
                    final_step_is_sulfonamide = True
                    findings_json["atomic_checks"]["named_reactions"].append("Sulfonamide synthesis (Schotten-Baumann) secondary amine")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
                            "position": "last_stage"
                        }
                    })
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children with adjusted depth
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage sulfonamide coupling detected: {final_step_is_sulfonamide}")
    return final_step_is_sulfonamide, findings_json
