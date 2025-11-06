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


# Refactored list of deprotection reactions
LATE_STAGE_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Alcohol deprotection from silyl ethers",
    "Alcohol deprotection from silyl ethers (double)",
    "Alcohol deprotection from silyl ethers (diol)",
    "Acetal hydrolysis to diol",
    "Acetal hydrolysis to aldehyde",
    "Ketal hydrolysis to ketone",
    "Hydroxyl benzyl deprotection",
    "Carboxyl benzyl deprotection",
    "Cleavage of methoxy ethers to alcohols",
    "Cleavage of alkoxy ethers to alcohols",
    "Ether cleavage to primary alcohol",
    "N-glutarimide deprotection",
    "Phthalimide deprotection",
    "TMS deprotection from alkyne",
    "Tert-butyl deprotection of amine",
    "COOH ethyl deprotection",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "Deprotection of carboxylic acid",
    "Hydrolysis of amides/imides/carbamates",
    "Hydrogenolysis of amides/imides/carbamates",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks if the final step (depth 1) of a synthesis is one of a specific set of deprotection reactions.
    The reactions checked for are defined in the LATE_STAGE_DEPROTECTION_REACTIONS list and include various Boc, silyl, acetal/ketal, benzyl, and ester deprotections.
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

    final_step_is_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_deprotection, findings_json

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # The final reaction step is at depth 1 (one step before the final product)
        if node["type"] == "reaction" and depth == 1:
            print(f"Analyzing potential final reaction at depth {depth}")

            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Found reaction SMILES: {rsmi}")

                # Check for common deprotection reactions using checker functions
                for reaction_type in LATE_STAGE_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        final_step_is_deprotection = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        print(f"Final step is {reaction_type}")
                        break

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {final_step_is_deprotection}")

    # Add structural constraint if the condition is met
    if final_step_is_deprotection:
        # This corresponds to the 'positional' constraint in the strategy JSON
        # We need to reconstruct the constraint object as it was in the original strategy JSON
        # The target list is the same as LATE_STAGE_DEPROTECTION_REACTIONS
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": LATE_STAGE_DEPROTECTION_REACTIONS,
                "position": "last_stage"
            }
        })

    return final_step_is_deprotection, findings_json
