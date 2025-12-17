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


DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Hydroxyl benzyl deprotection",
    "Carboxyl benzyl deprotection",
    "Cleavage of methoxy ethers to alcohols",
    "Cleavage of alkoxy ethers to alcohols",
    "Ether cleavage to primary alcohol",
    "COOH ethyl deprotection",
    "Tert-butyl deprotection of amine",
    "TMS deprotection from alkyne",
    "N-glutarimide deprotection",
    "Phthalimide deprotection",
    "Alcohol deprotection from silyl ethers",
    "Alcohol deprotection from silyl ethers (double)",
    "Alcohol deprotection from silyl ethers (diol)",
    "Acetal hydrolysis to diol",
    "Acetal hydrolysis to aldehyde",
    "Ketal hydrolysis to ketone",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the final synthetic step is a deprotection by matching the reaction against a curated list of named deprotection reactions specified in DEPROTECTION_REACTIONS.
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

    # In a retrosynthetic route, the root is the target molecule
    # The final synthetic step is the first reaction node we encounter

    def find_final_synthetic_step(node, depth=0):
        # print(f"Examining node at depth {depth}: {node.get('type', 'unknown')}")

        # If this is a reaction node and we're at depth 0 or 1, it's the final synthetic step
        if node["type"] == "reaction" and depth <= 1:
            # print(f"Found final synthetic step at depth {depth}")
            return node

        # If this is a molecule node, check its children
        if node["type"] == "mol" and "children" in node:
            for child in node["children"]:
                # Determine the new depth based on the current node's type
                new_depth = depth + 1 if node["type"] != "reaction" else depth
                result = find_final_synthetic_step(child, new_depth)
                if result:
                    return result

        return None

    # Find the final synthetic step
    final_step = find_final_synthetic_step(route)

    if not final_step:
        # print("Could not identify the final synthetic step")
        return result, findings_json

    # Check if the reaction is a deprotection
    if "metadata" in final_step and "rsmi" in final_step["metadata"]:
        rsmi = final_step["metadata"]["rsmi"]
        # print(f"Checking if final step is a deprotection: {rsmi}")

        # Check if the reaction matches any known deprotection reaction
        for reaction_type in DEPROTECTION_REACTIONS:
            if checker.check_reaction(reaction_type, rsmi):
                # print(f"Detected late-stage deprotection: {reaction_type}")
                result = True
                findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "deprotection_reaction",
                        "position": "last_stage"
                    }
                })
                break

    # if not result:
        # print("Final step is not a deprotection")
    return result, findings_json
