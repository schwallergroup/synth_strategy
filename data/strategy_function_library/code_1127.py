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


AMIDE_FORMATION_REACTION_TYPES = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis employs a late-stage amide formation strategy. It identifies an amide formation step by checking if the reaction matches a predefined list of named reactions (e.g., 'Carboxylic acid with primary amine to amide', 'Schotten-Baumann_amide'). A formation is considered 'late-stage' if it occurs within the final third of the synthetic route.
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

    amide_formation_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_depth, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                is_amide_formation = False
                for reaction_type in AMIDE_FORMATION_REACTION_TYPES:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_amide_formation = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                if is_amide_formation:
                    # In retrosynthesis, lower depth means later stage in forward synthesis
                    if amide_formation_depth is None or depth < amide_formation_depth:
                        amide_formation_depth = depth

            except Exception as e:
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reactions)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # In retrosynthesis, late stage means close to the target (low depth)
    # Define late stage as occurring in the first third of the max possible depth
    if max_depth == 0: # Handle single-step synthesis
        result = amide_formation_depth is not None
        if result:
            # Add the structural constraint if it's a single-step amide formation
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target_group": AMIDE_FORMATION_REACTION_TYPES,
                    "condition": "any",
                    "position": "within_final_third",
                    "description": "An amide formation reaction must occur within the final third of the synthesis (retrosynthesis depth <= max_depth / 3)."
                }
            })
        return result, findings_json

    is_late_stage = amide_formation_depth is not None and amide_formation_depth <= max_depth / 3

    if is_late_stage:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target_group": AMIDE_FORMATION_REACTION_TYPES,
                "condition": "any",
                "position": "within_final_third",
                "description": "An amide formation reaction must occur within the final third of the synthesis (retrosynthesis depth <= max_depth / 3)."
            }
        })

    return is_late_stage, findings_json
