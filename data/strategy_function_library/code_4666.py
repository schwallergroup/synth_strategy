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


BOC_PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
]

BOC_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a Boc protection/deprotection strategy by identifying reactions from predefined lists of named Boc protection and deprotection reaction types. A valid strategy is flagged if a protection reaction occurs earlier in the synthesis (higher depth value) than a subsequent deprotection reaction.
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

    # Track protection and deprotection events
    protection_events = []  # (molecule_smiles, depth)
    deprotection_events = []  # (molecule_smiles, depth)

    def dfs_traverse(node, depth=0):
        nonlocal findings_json, protection_events, deprotection_events
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            # Check for Boc protection reactions by name
            for rxn_type in BOC_PROTECTION_REACTIONS:
                if checker.check_reaction(rxn_type, rsmi):
                    protection_events.append((product_str, depth))
                    if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    break

            # Check for Boc deprotection reactions by name
            for rxn_type in BOC_DEPROTECTION_REACTIONS:
                if checker.check_reaction(rxn_type, rsmi):
                    reactants = reactants_str.split('.')
                    # Find which reactant had the Boc group
                    for reactant in reactants:
                        if checker.check_fg("Boc", reactant):
                            deprotection_events.append((reactant, depth))
                            if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Boc")
                            break # Assume only one reactant is deprotected
                    if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    break

        # Recursively process children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    result = False

    # Check if we have both protection and deprotection events
    if not protection_events or not deprotection_events:
        return result, findings_json

    # In retrosynthesis, protection should be at higher depth than deprotection
    # Find any valid protection-deprotection pair
    for _, protection_depth in protection_events:
        for _, deprotection_depth in deprotection_events:
            if protection_depth > deprotection_depth:
                result = True
                # Add the structural constraint if the condition is met
                structural_constraint_obj = {
                    "type": "sequence",
                    "details": {
                        "before": {
                            "type": "named_reaction",
                            "any_of": [
                                "Boc amine protection",
                                "Boc amine protection explicit",
                                "Boc amine protection with Boc anhydride",
                                "Boc amine protection (ethyl Boc)",
                                "Boc amine protection of secondary amine",
                                "Boc amine protection of primary amine"
                            ]
                        },
                        "after": {
                            "type": "named_reaction",
                            "any_of": [
                                "Boc amine deprotection",
                                "Boc amine deprotection of guanidine",
                                "Boc amine deprotection to NH-NH2",
                                "Tert-butyl deprotection of amine"
                            ]
                        }
                    }
                }
                if structural_constraint_obj not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(structural_constraint_obj)
                # No break here, as we want to find all valid pairs, but the result is already True

    return result, findings_json
