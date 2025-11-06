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


LATE_STAGE_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "COOH ethyl deprotection",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "Hydroxyl benzyl deprotection",
    "Carboxyl benzyl deprotection",
    "Cleavage of methoxy ethers to alcohols",
    "Cleavage of alkoxy ethers to alcohols",
    "Ether cleavage to primary alcohol",
    "N-glutarimide deprotection",
    "Phthalimide deprotection",
    "TMS deprotection from alkyne",
    "Tert-butyl deprotection of amine",
    "Deprotection of carboxylic acid",
    "Alcohol deprotection from silyl ethers",
    "Alcohol deprotection from silyl ethers (double)",
    "Alcohol deprotection from silyl ethers (diol)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for specific late-stage deprotection reactions at the final step of the synthesis. The reactions checked are: 'Boc amine deprotection', 'Boc amine deprotection of guanidine', 'Boc amine deprotection to NH-NH2', 'COOH ethyl deprotection', 'Ester saponification (methyl deprotection)', 'Ester saponification (alkyl deprotection)', 'Hydroxyl benzyl deprotection', 'Carboxyl benzyl deprotection', 'Cleavage of methoxy ethers to alcohols', 'Cleavage of alkoxy ethers to alcohols', 'Ether cleavage to primary alcohol', 'N-glutarimide deprotection', 'Phthalimide deprotection', 'TMS deprotection from alkyne', 'Tert-butyl deprotection of amine', 'Deprotection of carboxylic acid', 'Alcohol deprotection from silyl ethers', 'Alcohol deprotection from silyl ethers (double)', 'Alcohol deprotection from silyl ethers (diol)'.
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

    has_late_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_deprotection, findings_json

        if node["type"] == "reaction" and depth <= 1:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is a deprotection reaction
                for reaction_type in LATE_STAGE_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        has_late_deprotection = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        # Add the structural constraint if a deprotection reaction is found at the final/penultimate stage
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "any_listed_deprotection_reaction",
                                "position": "last_or_penultimate_stage"
                            }
                        })
                        return

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is 'chemical'
            next_depth = depth + 1

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    dfs_traverse(route)
    return has_late_deprotection, findings_json
