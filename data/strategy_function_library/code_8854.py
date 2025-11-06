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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves a late-stage esterification of a carboxylic acid.
    Late stage means at depth 0 or 1 (final or penultimate step).
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

    esterification_found = False
    esterification_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal esterification_found, esterification_depth, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for specific, valid esterification reactions of carboxylic acids.
                # Transesterification is excluded as it does not start from a carboxylic acid.
                # Hydrolysis/saponification is excluded as it's the reverse reaction.
                is_esterification = checker.check_reaction(
                    "Esterification of Carboxylic Acids", rsmi
                )
                is_o_alkylation = checker.check_reaction(
                    "O-alkylation of carboxylic acids with diazo compounds", rsmi
                )

                if is_esterification:
                    esterification_found = True
                    esterification_depth = depth
                    if "Esterification of Carboxylic Acids" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Esterification of Carboxylic Acids")
                
                if is_o_alkylation:
                    esterification_found = True
                    esterification_depth = depth
                    if "O-alkylation of carboxylic acids with diazo compounds" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("O-alkylation of carboxylic acids with diazo compounds")

            except Exception:
                # Ignore reactions that fail to parse
                pass

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is 'chemical', depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if esterification was found at depth 0 or 1 (late stage)
    result = (
        esterification_found and (esterification_depth is not None) and (esterification_depth <= 1)
    )
    
    if result:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "targets": [
                    "Esterification of Carboxylic Acids",
                    "O-alkylation of carboxylic acids with diazo compounds"
                ],
                "position_type": "depth_from_root",
                "operator": "<=",
                "value": 1
            }
        })

    return result, findings_json
