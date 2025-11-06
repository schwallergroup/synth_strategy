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


# List of alcohol to halide conversion reaction types
ALCOHOL_TO_HALIDE_REACTIONS = [
    "Alcohol to chloride_sulfonyl chloride",
    "Alcohol to chloride_SOCl2",
    "Alcohol to chloride_PCl5_ortho",
    "Alcohol to chloride_POCl3_ortho",
    "Alcohol to chloride_POCl3_para",
    "Alcohol to chloride_POCl3",
    "Alcohol to chloride_HCl",
    "Appel reaction",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis includes a late-stage alcohol to halide conversion.
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

    found_alcohol_to_halide = False

    def dfs_traverse(node, depth=0):
        nonlocal found_alcohol_to_halide, findings_json

        if found_alcohol_to_halide:
            return

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check if this is a late-stage reaction (depth ≤ 2)
            if depth <= 2:
                # Check if this is a known alcohol-to-halide reaction type
                for reaction_type in ALCOHOL_TO_HALIDE_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        found_alcohol_to_halide = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "alcohol_to_halide_conversion",
                                "depth_constraint": {
                                    "operator": "<=",
                                    "value": 2
                                }
                            }
                        })
                        return

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return found_alcohol_to_halide, findings_json
