from typing import Tuple, Dict, List
import copy
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


# Refactored lists for enumeration
PROTECTION_REACTIONS = [
    "Alcohol protection with silyl ethers",
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
    "Protection of carboxylic acid",
    "Aldehyde or ketone acetalization",
    "Diol acetalization",
]

DEPROTECTION_REACTIONS = [
    "Alcohol deprotection from silyl ethers",
    "Alcohol deprotection from silyl ethers (double)",
    "Alcohol deprotection from silyl ethers (diol)",
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Deprotection of carboxylic acid",
    "Acetal hydrolysis to diol",
    "Acetal hydrolysis to aldehyde",
    "Ketal hydrolysis to ketone",
    "Tert-butyl deprotection of amine",
    "Hydroxyl benzyl deprotection",
    "Carboxyl benzyl deprotection",
    "Cleavage of methoxy ethers to alcohols",
    "Cleavage of alkoxy ethers to alcohols",
    "Ether cleavage to primary alcohol",
    "COOH ethyl deprotection",
    "N-glutarimide deprotection",
    "Phthalimide deprotection",
    "TMS deprotection from alkyne",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis involves both protection and deprotection steps.
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

    protection_found = False
    deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal protection_found, deprotection_found, findings_json

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for protection reactions if not already found
                if not protection_found:
                    for reaction_type in PROTECTION_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            protection_found = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            break

                # Check for deprotection reactions if not already found
                if not deprotection_found:
                    for reaction_type in DEPROTEPTION_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            deprotection_found = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            break

        # Determine the new depth for recursive calls
        new_depth = depth
        if node["type"] != "reaction":  # If current node is 'chemical'
            new_depth = depth + 1

        # Continue traversing the route
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)
            # If we've already found both, we can stop traversing
            if protection_found and deprotection_found:
                return

    dfs_traverse(route, 0)
    
    result = protection_found and deprotection_found
    
    if result:
        # Add the structural constraint if both protection and deprotection steps are found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "protection_reaction",
                    "deprotection_reaction"
                ]
            }
        })

    return result, findings_json
