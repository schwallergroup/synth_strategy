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


DEPROTECTION_REACTION_TYPES = [
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
]

PROTECTING_GROUPS_OF_INTEREST = [
    "Boc",
    "TMS ether protective group",
    "Silyl protective group",
    "Acetal/Ketal",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage deprotection strategy. This is identified if the final synthetic step (depth=1) matches a reaction in the `DEPROTECTION_REACTION_TYPES` list, or if it involves the removal of a functional group specified in the `PROTECTING_GROUPS_OF_INTEREST` list.
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

    has_final_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_final_deprotection, findings_json

        # Check if this is the final reaction (depth 1 leads to final product at depth 0)
        if node["type"] == "reaction" and depth == 1:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check if this is any type of deprotection reaction
            deprotection_found_by_reaction_type = False
            for reaction_type in DEPROTECTION_REACTION_TYPES:
                if checker.check_reaction(reaction_type, rsmi):
                    has_final_deprotection = True
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    deprotection_found_by_reaction_type = True
                    break
            
            if deprotection_found_by_reaction_type:
                # Add the structural constraint if a deprotection reaction was found at depth 1
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "deprotection_reaction",
                        "position": "last_stage"
                    }
                })

            # If no specific reaction type matched, check for general deprotection pattern
            if not has_final_deprotection:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for various protection groups
                for fg in PROTECTING_GROUPS_OF_INTEREST:
                    fg_in_reactants = any(checker.check_fg(fg, r) for r in reactants_smiles)
                    fg_in_product = checker.check_fg(fg, product_smiles)

                    if fg_in_reactants and not fg_in_product:
                        has_final_deprotection = True
                        findings_json["atomic_checks"]["functional_groups"].append(fg)
                        findings_json["atomic_checks"]["named_reactions"].append("protecting_group_removal") # As per strategy JSON
                        # Add the structural constraint if a deprotection by FG removal was found at depth 1
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "deprotection_reaction",
                                "position": "last_stage"
                            }
                        })
                        break

        # Process children
        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return has_final_deprotection, findings_json
