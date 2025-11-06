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


NITRO_ADDITION_REACTIONS_RETRO = [
    "Reduction of nitro groups to amines",
    "Nef reaction (nitro to ketone)",
    "Aromatic nitration with HNO3",
    "Aromatic nitration with NO3 salt",
    "Aromatic nitration with NO2 salt",
    "Aromatic nitration with alkyl NO2",
    "Non-aromatic nitration with HNO3",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis involves late-stage nitro group removal (at depth 0 or 1).
    In retrosynthesis, this means looking for nitro group addition reactions at depths 0 or 1.
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

    has_late_stage_nitro_removal = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_nitro_removal, findings_json

        if node["type"] == "reaction" and depth <= 1:
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    product_smiles = rsmi.split(">")[-1]

                    if checker.check_fg("Nitro group", product_smiles):
                        if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

                        for reaction_type in NITRO_ADDITION_REACTIONS_RETRO:
                            if checker.check_reaction(reaction_type, rsmi):
                                if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                                has_late_stage_nitro_removal = True
                                # Add structural constraint if both conditions are met
                                if {"type": "positional", "details": {"target": "A nitro addition reaction where the product contains a nitro group", "position": "last or penultimate stage"}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "A nitro addition reaction where the product contains a nitro group", "position": "last or penultimate stage"}})
                                return
            except Exception:
                pass

        if not has_late_stage_nitro_removal:
            for child in node.get("children", []):
                # Determine the new depth based on the current node's type
                new_depth = depth
                if node["type"] != "reaction":  # If current node is chemical, depth increases
                    new_depth = depth + 1
                
                dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return has_late_stage_nitro_removal, findings_json
