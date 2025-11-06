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


BOC_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if Boc protection is maintained throughout the synthesis.
    Returns True if Boc protection is added early in the synthesis and
    maintained until the final product without deprotection steps.
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

    # Track Boc protection status and deprotection events
    boc_status = {
        "has_boc_at_depth": {},  # Track depths where Boc is present
        "deprotection_occurred": False,  # Flag if Boc deprotection occurred
        "max_depth": 0,  # Maximum depth reached
        "boc_in_final": False,  # If final product has Boc
    }

    def dfs_traverse(node, depth=0):
        nonlocal boc_status, findings_json

        boc_status["max_depth"] = max(boc_status["max_depth"], depth)

        if node["type"] == "mol":
            # Check if molecule has Boc group
            if checker.check_fg("Boc", node["smiles"]):
                boc_status["has_boc_at_depth"][depth] = True
                if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boc")
                print(f"Found Boc protection at depth {depth}")

                # If this is the final product (depth 0), mark it
                if depth == 0:
                    boc_status["boc_in_final"] = True

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Check if this is a Boc deprotection reaction
            rxn_smiles = node["metadata"]["mapped_reaction_smiles"]
            for name in BOC_DEPROTECTION_REACTIONS:
                if checker.check_reaction(name, rxn_smiles):
                    print(f"Found Boc deprotection reaction at depth {depth}")
                    boc_status["deprotection_occurred"] = True
                    if name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                    break

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'mol' (chemical)
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if Boc is maintained throughout
    # Conditions:
    # 1. Boc must be present in the final product
    # 2. No Boc deprotection should have occurred
    # 3. Boc should be present at some early stage (depth > 0)

    early_stage_boc = any(
        boc_status["has_boc_at_depth"].get(d, False) for d in range(1, boc_status["max_depth"] + 1)
    )

    maintained = (
        boc_status["boc_in_final"] and early_stage_boc and not boc_status["deprotection_occurred"]
    )

    print(f"Boc in final: {boc_status['boc_in_final']}")
    print(f"Early stage Boc: {early_stage_boc}")
    print(f"Deprotection occurred: {boc_status['deprotection_occurred']}")
    print(f"Boc maintained: {maintained}")

    # Record structural constraints if conditions are met
    if boc_status["boc_in_final"]:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Boc",
                "target_type": "functional_group",
                "position": "last_stage"
            }
        })
    
    if early_stage_boc:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Boc",
                "target_type": "functional_group",
                "position": "not_last_stage"
            }
        })

    if not boc_status["deprotection_occurred"]:
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "targets": [
                    "Boc amine deprotection",
                    "Boc amine deprotection of guanidine",
                    "Boc amine deprotection to NH-NH2",
                    "Tert-butyl deprotection of amine"
                ],
                "target_type": "named_reaction"
            }
        })

    return maintained, findings_json
