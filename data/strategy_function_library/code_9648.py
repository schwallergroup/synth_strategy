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
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Ester with ammonia to amide",
    "Acyl chloride with ammonia to amide",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis involves a late-stage amide bond formation by checking against a predefined list of common amide coupling reactions.
    Late-stage is defined as occurring at depth 0 or 1 in the synthesis route.
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

    amide_formation_at_depth = {}

    def dfs_traverse(node, depth=0):
        nonlocal findings_json, amide_formation_at_depth
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check if this is a known amide formation reaction
            is_amide_formation = False
            for reaction_name in AMIDE_FORMATION_REACTIONS:
                if checker.check_reaction(reaction_name, rsmi):
                    is_amide_formation = True
                    # Record the detected reaction in findings_json
                    if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    break

            # Record amide formation if it's a known amide formation reaction
            if is_amide_formation:
                amide_formation_at_depth[depth] = True
                print(f"Amide formation detected at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                print("Recognized as amide formation reaction")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Check if amide formation occurs at depth 0 or 1 (late stage)
    result = any(depth <= 1 for depth in amide_formation_at_depth.keys())

    if result:
        # Add the structural constraint if late-stage amide formation is detected
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "amide_formation",
                "position": "late_stage (depth <= 1)"
            }
        })

    print(f"Late-stage amide formation detected: {result}")
    print(f"Amide formations by depth: {amide_formation_at_depth}")

    return result, findings_json
