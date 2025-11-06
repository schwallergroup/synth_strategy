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


AMIDE_FORMATION_REACTIONS = [
    "Carboxylic acid with primary amine to amide",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects late-stage amide formation (depth ≤ 2).
    It looks for a reaction where a carboxylic acid and an amine form an amide bond.
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

    late_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_amide_formation, findings_json

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check if this is a late-stage reaction (depth ≤ 2)
            if depth <= 2:
                # Check for amide formation reactions directly
                for reaction_type in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        late_amide_formation = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        # Add the structural constraint if the condition is met
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "amide_formation",
                                "position": "late_stage (depth <= 2)"
                            }
                        })
                        return

        for child in node.get("children", []):
            if (
                not late_amide_formation
            ):  # Stop traversal if we already found what we're looking for
                new_depth = depth
                if node["type"] != "reaction": # If current node is chemical, depth increases
                    new_depth = depth + 1
                # If current node is reaction, depth remains the same
                dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return late_amide_formation, findings_json
