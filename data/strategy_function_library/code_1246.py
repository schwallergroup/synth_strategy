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


PROTECTION_REACTIONS_OF_INTEREST = [
    "Alcohol protection with silyl ethers",
    "Protection of carboxylic acid",
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
    "Aldehyde or ketone acetalization",
    "Diol acetalization",
]

DEPROTECTION_REACTIONS_OF_INTEREST = [
    "Acetal hydrolysis to diol",
    "Acetal hydrolysis to aldehyde",
    "Ketal hydrolysis to ketone",
    "Alcohol deprotection from silyl ethers",
    "Alcohol deprotection from silyl ethers (double)",
    "Alcohol deprotection from silyl ethers (diol)",
    "Deprotection of carboxylic acid",
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "TMS deprotection from alkyne",
    "Tert-butyl deprotection of amine",
    "Phthalimide deprotection",
    "N-glutarimide deprotection",
    "Hydroxyl benzyl deprotection",
    "Carboxyl benzyl deprotection",
    "Cleavage of methoxy ethers to alcohols",
    "Cleavage of alkoxy ethers to alcohols",
    "Ether cleavage to primary alcohol",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves specific, enumerated protection and deprotection reactions by checking against predefined lists of known transformations for common protecting groups like silyl ethers, Boc, and acetals.
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

    protection_reactions = []
    deprotection_reactions = []

    def find_protection_deprotection(node, depth=0):
        nonlocal protection_reactions, deprotection_reactions, findings_json
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rxn_smiles = node["metadata"]["mapped_reaction_smiles"]

            # Check for protection reactions
            for r in PROTECTION_REACTIONS_OF_INTEREST:
                if checker.check_reaction(r, rxn_smiles):
                    protection_reactions.append((rxn_smiles, depth))
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)

            # Check for deprotection reactions
            for r in DEPROTECTION_REACTIONS_OF_INTEREST:
                if checker.check_reaction(r, rxn_smiles):
                    deprotection_reactions.append((rxn_smiles, depth))
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            find_protection_deprotection(child, new_depth)

    find_protection_deprotection(route)

    # Remove duplicates
    unique_protection = set(rxn for rxn, _ in protection_reactions)
    unique_deprotection = set(rxn for rxn, _ in deprotection_reactions)

    # Consider it a protection/deprotection strategy if at least one protection and one deprotection are found
    has_protection = len(unique_protection) > 0
    has_deprotection = len(unique_deprotection) > 0
    print(
        f"Protection reactions: {len(unique_protection)}, Deprotection reactions: {len(unique_deprotection)}"
    )

    result = has_protection or has_deprotection

    if result:
        # Add the structural constraint if any protection or deprotection reaction was found
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "any_protection_or_deprotection_reaction",
                "operator": ">=",
                "value": 1
            }
        })

    return result, findings_json