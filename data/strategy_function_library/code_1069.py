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


O_ALKYLATION_REACTIONS = [
    "Williamson Ether Synthesis",
    "Williamson Ether Synthesis (intra to epoxy)",
    "O-alkylation of carboxylic acids with diazo compounds",
    "O-alkylation of amides with diazo compounds",
    "Mitsunobu aryl ether",
    "Chan-Lam etherification",
    "O-methylation",
    "Methylation of OH with DMS",
]


def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a late-stage O-alkylation occurs, defined as a reaction within the final synthetic step (depth=1). It identifies this by checking for a predefined list of O-alkylation named reactions, as defined in the O_ALKYLATION_REACTIONS list.
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

    late_stage_o_alkylation_detected = False

    def dfs_traverse(node, current_depth=0):
        nonlocal late_stage_o_alkylation_detected, findings_json

        node["depth"] = current_depth

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            depth = node.get("depth", None)
            if depth <= 1:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                for reaction_name in O_ALKYLATION_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        late_stage_o_alkylation_detected = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        # Add the structural constraint if a late-stage O-alkylation is detected
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "targets": [
                                    "Williamson Ether Synthesis",
                                    "Williamson Ether Synthesis (intra to epoxy)",
                                    "O-alkylation of carboxylic acids with diazo compounds",
                                    "O-alkylation of amides with diazo compounds",
                                    "Mitsunobu aryl ether",
                                    "Chan-Lam etherification",
                                    "O-methylation",
                                    "Methylation of OH with DMS"
                                ],
                                "position": "last_stage"
                            }
                        })
                        break

        # Traverse children with incremented depth
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemical nodes)
                dfs_traverse(child, current_depth)
            else:
                # If current node is not a reaction (e.g., chemical), depth increases for children (reaction nodes)
                dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)
    return late_stage_o_alkylation_detected, findings_json
