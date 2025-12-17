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


ALCOHOL_TO_HALIDE_REACTIONS = [
    "Alcohol to chloride",
    "Alcohol to bromide",
    "Alcohol to iodide",
    "Alcohol to chloride_sulfonyl chloride",
    "Alcohol to chloride_SOCl2",
    "Alcohol to chloride_CHCl3",
    "Alcohol to chloride_CH2Cl2",
    "Alcohol to chloride_PCl5_ortho",
    "Alcohol to chloride_POCl3_ortho",
    "Alcohol to chloride_POCl3_para",
    "Alcohol to chloride_POCl3",
    "Alcohol to chloride_HCl",
    "Alcohol to chloride_Salt",
    "Alcohol to chloride_Other",
    "Appel reaction",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a hydroxyl group is converted to a halide early in the synthesis.

    This strategy is identified by checking if a reaction step matches a known
    alcohol-to-halide transformation from a predefined list (e.g., Appel reaction,
    reaction with SOCl2) and occurs in the first half of the synthesis sequence
    (i.e., depth > max_depth / 2).
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

    hydroxyl_to_halide_found = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal hydroxyl_to_halide_found, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    return

                # Check if this is a known alcohol to halide reaction
                is_alcohol_to_halide_rxn = False
                for rxn in ALCOHOL_TO_HALIDE_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        is_alcohol_to_halide_rxn = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break

                if is_alcohol_to_halide_rxn:
                    # Early in synthesis = higher depth in retrosynthetic tree
                    if depth > max_depth / 2:
                        hydroxyl_to_halide_found = True
                        # Add structural constraint if the condition is met
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "alcohol_to_halide_conversion",
                                "position": "first_half"
                            }
                        })
            except Exception:
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains same when going from reaction to chemical
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # First pass to determine max depth
    def calculate_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            # This max_depth calculation should follow the same new logic
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            calculate_max_depth(child, new_depth)

    calculate_max_depth(route)

    # Second pass to find the strategy
    dfs_traverse(route)

    return hydroxyl_to_halide_found, findings_json
