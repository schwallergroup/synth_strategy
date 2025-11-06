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
    This function detects if the synthesis involves a final oxidation step
    that converts an aldehyde to a carboxylic acid.
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

    final_oxidation = False

    def dfs_traverse(node, depth=0):
        nonlocal final_oxidation, findings_json

        # Only analyze the final reaction step (depth=1 in retrosynthesis traversal)
        if node["type"] == "reaction" and depth == 1:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                # The most robust check is to use the specific reaction template.
                # The original code included checks for other oxidation types and
                # brittle string searches for oxidants, which are sources of
                # false positives and have been removed.
                if checker.check_reaction(
                    "Oxidation of aldehydes to carboxylic acids", rsmi
                ):
                    final_oxidation = True
                    findings_json["atomic_checks"]["named_reactions"].append("Oxidation of aldehydes to carboxylic acids")
                    # Add the structural constraint if the condition is met
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Oxidation of aldehydes to carboxylic acids",
                            "position": "last_stage"
                        }
                    })
            except Exception:
                # Ignore potential errors in reaction analysis to avoid crashing.
                pass

        for child in node.get("children", []):
            # Optimization: stop traversing if we've found our match.
            # This is a valid control flow modification as it doesn't change the logic, only efficiency.
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            
            if final_oxidation:
                break
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return final_oxidation, findings_json
