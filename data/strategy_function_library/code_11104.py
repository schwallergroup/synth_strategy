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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

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


LATE_STAGE_AMINE_ALKYLATIONS = [
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for late-stage (depth <= 2) amine incorporation using a defined list of N-alkylation reactions. The specific reactions checked are: 'N-alkylation of primary amines with alkyl halides' and 'N-alkylation of secondary amines with alkyl halides'.
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

    amine_incorporation = False

    def dfs_traverse(node, depth=0):
        nonlocal amine_incorporation, findings_json

        if (
            node["type"] == "reaction"
            and depth <= 2
            and "metadata" in node
            and "mapped_reaction_smiles" in node["metadata"]
        ):
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check if the reaction is one of the specified N-alkylation types
            for rxn_name in LATE_STAGE_AMINE_ALKYLATIONS:
                if checker.check_reaction(rxn_name, rsmi):
                    amine_incorporation = True
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                    # Add structural constraint if depth condition is met
                    if depth <= 2:
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": [
                                    "N-alkylation of primary amines with alkyl halides",
                                    "N-alkylation of secondary amines with alkyl halides"
                                ],
                                "position": {
                                    "constraint": "depth",
                                    "operator": "<=",
                                    "value": 2
                                }
                            }
                        })
                    return

        for child in node.get("children", []):
            # Stop traversing if we've already found our match
            if amine_incorporation:
                break
            
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    return amine_incorporation, findings_json
