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


AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Acyl chloride with ammonia to amide",
    "Ester with ammonia to amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the final step (depth 1) is one of several common amide coupling reactions, as defined in the AMIDE_COUPLING_REACTIONS list.
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

    final_step_is_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_amide_coupling, findings_json

        # Check if this is a reaction node at the final step (depth 1)
        if node["type"] == "reaction" and depth == 1:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is an amide coupling using reaction checkers
                for reaction_type in AMIDE_COUPLING_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        final_step_is_amide_coupling = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        # Since we found a reaction at depth 1, this implies the structural constraint
                        # is met if final_step_is_amide_coupling remains True after traversal.
                        # We will add the structural constraint after the traversal.
                        return
            except Exception:
                # Error analyzing, so we cannot confirm the strategy for this node
                pass

        # Traverse children with incremented depth
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., chemical), increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # After traversal, if the flag is set, add the structural constraint
    if final_step_is_amide_coupling:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "amide_coupling",
                "position": "last_stage"
            }
        })

    return final_step_is_amide_coupling, findings_json
