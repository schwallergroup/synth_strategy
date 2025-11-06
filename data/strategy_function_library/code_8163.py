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
    This function detects if the synthesis involves a late-stage
    nitro reduction to form an amine.
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

    # Track if nitro reduction occurs at late stage (depth ≤ 2)
    late_stage_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_reduction, findings_json

        if node["type"] == "reaction" and depth <= 2:  # Late stage (depth 0-2)
            # Extract reactants and product
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                # Primary check: Use the reaction checker for accuracy.
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    late_stage_reduction = True
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                    # Add structural constraint if late_stage_reduction is true
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Reduction of nitro groups to amines",
                            "position": "late_stage"
                        }
                    })
                    return

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only from chemical to reaction node
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    print(f"Late stage nitro reduction detected: {late_stage_reduction}")
    return late_stage_reduction, findings_json
