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
    Detects a Boc-based protecting group strategy. A strategy is identified if either of two conditions is met: 1) A Boc group is present on at least one molecule in the route and no Boc deprotection reactions occur. 2) The route contains at least one Boc protection reaction and at least one Boc deprotection reaction. The specific reactions are defined in the module-level lists `BOC_PROTECTION_REACTIONS` and `BOC_DEPROTECTION_REACTIONS`.
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

    # Track if we've seen Boc protection/deprotection reactions
    boc_protection_reaction = False
    boc_deprotection_reaction = False

    # Track molecules with Boc groups
    molecules_with_boc = {}

    result = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_reaction, boc_deprotection_reaction, findings_json

        if node["type"] == "mol":
            # Check if this molecule has a Boc group
            if checker.check_fg("Boc", node["smiles"]):
                molecules_with_boc[node["smiles"]] = depth
                if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boc")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check if this is a Boc protection or deprotection reaction
            for r in BOC_PROTECTION_REACTIONS:
                if checker.check_reaction(r, rsmi):
                    boc_protection_reaction = True
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)

            for r in BOC_DEPROTECTION_REACTIONS:
                if checker.check_reaction(r, rsmi):
                    boc_deprotection_reaction = True
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children
                dfs_traverse(child, depth)
            else:
                # If current node is not a reaction (e.g., chemical), depth increases
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have molecules with Boc groups
    has_boc_molecules = len(molecules_with_boc) > 0

    # A true Boc protection strategy would have:
    # 1. Boc groups present in molecules
    # 2. Either:
    #    a. No deprotection reactions (Boc maintained throughout)
    #    b. Protection and deprotection reactions both occur in the route

    if has_boc_molecules:
        if not boc_deprotection_reaction:
            # Boc is present and never removed - clear protection strategy
            result = True
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "description": "A Boc functional group is present on at least one molecule, and no Boc deprotection reactions occur anywhere in the route.",
                    "present_fgs": [
                        "Boc"
                    ],
                    "absent_reactions": [
                        "Boc amine deprotection",
                        "Boc amine deprotection of guanidine",
                        "Boc amine deprotection to NH-NH2",
                        "Tert-butyl deprotection of amine"
                    ]
                }
            })
        elif boc_protection_reaction and boc_deprotection_reaction:
            # A full protect-deprotect sequence exists
            result = True
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "description": "The route contains at least one reaction from the 'Boc protection' group and at least one reaction from the 'Boc deprotection' group.",
                    "present_from_groups": [
                        {
                            "group_name": "Boc protection reactions",
                            "members": [
                                "Boc amine protection",
                                "Boc amine protection explicit",
                                "Boc amine protection with Boc anhydride",
                                "Boc amine protection (ethyl Boc)",
                                "Boc amine protection of secondary amine",
                                "Boc amine protection of primary amine"
                            ]
                        },
                        {
                            "group_name": "Boc deprotection reactions",
                            "members": [
                                "Boc amine deprotection",
                                "Boc amine deprotection of guanidine",
                                "Boc amine deprotection to NH-NH2",
                                "Tert-butyl deprotection of amine"
                            ]
                        }
                    ]
                }
            })

    return result, findings_json
