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


LATE_STAGE_N_MODIFICATIONS = [
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Methylation with MeI_primary",
    "Methylation with MeI_secondary",
    "N-methylation",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of secondary amines",
    "Acylation of primary amines",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the final synthetic step (depth=1) is a specific N-alkylation or N-acylation reaction. The check is performed by matching against a predefined list of reaction names, including various alkylations and acylations of primary and secondary amines.
    """
    found = False
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal found, findings_json

        # Process reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Only check the final step of the synthesis
            if depth == 1:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                
                # Check if the reaction matches any of the specified N-modification types
                for name in LATE_STAGE_N_MODIFICATIONS:
                    if checker.check_reaction(name, rsmi):
                        found = True
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                        # Add the structural constraint if found at the last stage
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "targets": [
                                    "N-alkylation of primary amines with alkyl halides",
                                    "N-alkylation of secondary amines with alkyl halides",
                                    "Methylation with MeI_primary",
                                    "Methylation with MeI_secondary",
                                    "N-methylation",
                                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                    "Acylation of secondary amines",
                                    "Acylation of primary amines"
                                ],
                                "position": "last_stage"
                            }
                        })
                        return # Strategy found, no need to check children of this node

        # Traverse children (depth-first)
        # Stop traversing if the condition has already been met
        if not found:
            for child in node.get("children", []):
                # New logic for depth calculation
                new_depth = depth
                if node["type"] != "reaction": # If current node is chemical, depth increases
                    new_depth = depth + 1
                # If current node is reaction, depth remains the same
                dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return found, findings_json
