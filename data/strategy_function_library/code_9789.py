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
    Detects a Boc protection/deprotection sequence by identifying specific named reactions. Protection reactions checked are 'Boc amine protection', 'Boc amine protection explicit', 'Boc amine protection with Boc anhydride', 'Boc amine protection (ethyl Boc)', 'Boc amine protection of secondary amine', and 'Boc amine protection of primary amine'. Deprotection reactions checked are 'Boc amine deprotection', 'Boc amine deprotection of guanidine', 'Boc amine deprotection to NH-NH2', and 'Tert-butyl deprotection of amine'.
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

    # Track all protection and deprotection steps
    protection_depths = []
    deprotection_depths = []
    
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal protection_depths, deprotection_depths, findings_json, result

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for Boc protection reactions
            for name in BOC_PROTECTION_REACTIONS:
                if checker.check_reaction(name, rsmi):
                    protection_depths.append(depth)
                    if name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                    break # Only need to find one match per reaction node

            # Check for Boc deprotection reactions
            for name in BOC_DEPROTECTION_REACTIONS:
                if checker.check_reaction(name, rsmi):
                    deprotection_depths.append(depth)
                    if name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                    break # Only need to find one match per reaction node

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if we have both protection and deprotection in the correct order
    # In forward synthesis: protection happens first (higher depth in retrosynthesis)
    # then deprotection happens later (lower depth in retrosynthesis)
    for p_depth in protection_depths:
        for d_depth in deprotection_depths:
            # In retrosynthetic traversal, protection should be at a higher depth than deprotection
            if p_depth > d_depth:
                result = True
                # Add the structural constraint if the condition is met
                structural_constraint_obj = {
                    "type": "sequence",
                    "details": {
                        "before": "Any Boc protection reaction from the checked list",
                        "after": "Any Boc deprotection reaction from the checked list",
                        "notes": "In the retrosynthesis tree, the protection step must occur at a greater depth than the deprotection step."
                    }
                }
                if structural_constraint_obj not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(structural_constraint_obj)
                # No need to break here, as we want to find all pairs, but the result is already True

    return result, findings_json