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


LATE_STAGE_ALCOHOL_TO_HALIDE_REACTIONS = [
    "Alcohol to chloride_SOCl2",
    "Alcohol to chloride_PCl5_ortho",
    "Alcohol to chloride_POCl3_ortho",
    "Alcohol to chloride_POCl3_para",
    "Alcohol to chloride_POCl3",
    "Alcohol to chloride_HCl",
    "Alcohol to chloride_Salt",
    "Alcohol to chloride_Other",
    "Alcohol to chloride_sulfonyl chloride",
    "Appel reaction",
    "Alcohol to bromide",
    "Alcohol to iodide",
    "Alcohol to fluoride",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks if the final synthetic step is a known alcohol-to-halide conversion. This is verified by checking the reaction
    against a curated list of relevant named reactions, including the Appel reaction and transformations using reagents
    like SOCl2, PCl5, and PBr3. The list of checked reactions is defined in `LATE_STAGE_ALCOHOL_TO_HALIDE_REACTIONS`.
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

    found_halogenation_at_final_step = False

    def dfs_traverse(node, depth=0):
        nonlocal found_halogenation_at_final_step, findings_json

        # Terminate early if we've already found the pattern
        if found_halogenation_at_final_step:
            return

        if node["type"] == "reaction" and depth == 1:
            try:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    return

                # Check if this is a known alcohol-to-halide reaction
                for reaction_type in LATE_STAGE_ALCOHOL_TO_HALIDE_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        found_halogenation_at_final_step = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        # Add the structural constraint if this is the final step and a relevant reaction is found
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "targets": [
                                    "Alcohol to chloride_SOCl2",
                                    "Alcohol to chloride_PCl5_ortho",
                                    "Alcohol to chloride_POCl3_ortho",
                                    "Alcohol to chloride_POCl3_para",
                                    "Alcohol to chloride_POCl3",
                                    "Alcohol to chloride_HCl",
                                    "Alcohol to chloride_Salt",
                                    "Alcohol to chloride_Other",
                                    "Alcohol to chloride_sulfonyl chloride",
                                    "Appel reaction",
                                    "Alcohol to bromide",
                                    "Alcohol to iodide",
                                    "Alcohol to fluoride"
                                ],
                                "position": "last_stage"
                            }
                        })
                        return  # Strategy found, no need to check further

            except Exception as e:
                # In case of an error, we just skip this node
                print(f"Error analyzing reaction: {e}")

        # Traverse children with new depth logic
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return found_halogenation_at_final_step, findings_json
