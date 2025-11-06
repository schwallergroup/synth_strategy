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


NITRATION_REACTION_TYPES = [
    "Aromatic nitration with HNO3",
    "Aromatic nitration with NO3 salt",
    "Aromatic nitration with NO2 salt",
    "Aromatic nitration with alkyl NO2",
    "Non-aromatic nitration with HNO3",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a late-stage functionalization strategy involving specific nitration reactions. It identifies reactions that match a predefined list of types (e.g., 'Aromatic nitration with HNO3') and confirms the net addition of a nitro group. The strategy is flagged only if it occurs in the final step of the synthesis (depth <= 1).
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

    nitration_detected = False
    nitration_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal nitration_detected, nitration_depth, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                is_nitration_reaction = False
                for r_type in NITRATION_REACTION_TYPES:
                    if checker.check_reaction(r_type, rsmi):
                        is_nitration_reaction = True
                        if r_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_type)

                if is_nitration_reaction:
                    has_nitro_in_product = checker.check_fg("Nitro group", product)
                    if has_nitro_in_product and "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

                    has_nitro_in_reactants = any(
                        checker.check_fg("Nitro group", r) for r in reactants
                    )

                    if has_nitro_in_product and not has_nitro_in_reactants:
                        nitration_detected = True
                        nitration_depth = min(nitration_depth, depth)

            except KeyError:
                pass
            except Exception:
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = nitration_detected and nitration_depth <= 1

    if result:
        # Add structural constraints if the overall strategy is detected
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "A nitration reaction from the predefined list that adds a nitro group",
                "position": "depth <= 1"
            }
        })
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "Nitro group in any reactant of the qualifying nitration step"
            }
        })

    return result, findings_json
