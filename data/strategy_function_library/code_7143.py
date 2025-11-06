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
rng_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**rng_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)

def dfs_traverse(node, depth):
    node["depth"] = depth
    if "children" in node:
        for child in node["children"]:
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same for its children
            dfs_traverse(child, new_depth)

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where fluorination is the final step
    in the synthesis route.
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
    result = False

    # The final product must have a reaction leading to it.
    if "children" not in route or not route["children"]:
        return result, findings_json

    # The node leading to the final product must be a reaction.
    last_reaction = route["children"][0]
    if (
        last_reaction["type"] != "reaction"
        or "metadata" not in last_reaction
        or "rsmi" not in last_reaction["metadata"]
    ):
        return result, findings_json

    reaction_smiles = last_reaction["metadata"]["rsmi"]

    # Rely on curated checker functions to identify fluorination reactions.
    is_aromatic_fluorination = checker.check_reaction(
        "Aromatic fluorination", reaction_smiles
    )
    is_fluorination = checker.check_reaction("Fluorination", reaction_smiles)

    if is_aromatic_fluorination:
        findings_json["atomic_checks"]["named_reactions"].append("Aromatic fluorination")
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Aromatic fluorination",
                "position": "last_stage"
            }
        })
        result = True

    if is_fluorination:
        findings_json["atomic_checks"]["named_reactions"].append("Fluorination")
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Fluorination",
                "position": "last_stage"
            }
        })
        result = True

    return result, findings_json
