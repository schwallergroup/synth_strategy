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
    Detects late-stage deprotection of a Boc-protected amine by checking for a specific list of named reactions. The reactions checked are: 'Boc amine deprotection', 'Boc amine deprotection of guanidine', 'Boc amine deprotection to NH-NH2', and 'Tert-butyl deprotection of amine'.
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

    boc_deprotection_at_late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_deprotection_at_late_stage, findings_json

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            print(f"Checking reaction at depth {depth}, RSMI: {rsmi}")

            is_boc_deprotection_found = False
            for rxn in BOC_DEPROTECTION_REACTIONS:
                if checker.check_reaction(rxn, rsmi):
                    is_boc_deprotection_found = True
                    if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                    break

            print(f"Is Boc deprotection: {is_boc_deprotection_found}")

            if is_boc_deprotection_found:
                # Consider depths 1, 2, and 3 as late-stage (root product is depth 0).
                if depth <= 3:
                    boc_deprotection_at_late_stage = True
                    print(f"Late-stage Boc deprotection detected at depth {depth}")
                    # Add the structural constraint if it's a late-stage deprotection
                    # This corresponds to the 'positional' constraint in the strategy JSON
                    if {"type": "positional", "details": {"targets": ["Boc amine deprotection", "Boc amine deprotection of guanidine", "Boc amine deprotection to NH-NH2", "Tert-butyl deprotection of amine"], "position": "late_stage", "definition": "depth <= 3"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "targets": [
                                    "Boc amine deprotection",
                                    "Boc amine deprotection of guanidine",
                                    "Boc amine deprotection to NH-NH2",
                                    "Tert-butyl deprotection of amine"
                                ],
                                "position": "late_stage",
                                "definition": "depth <= 3"
                            }
                        })

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node['type'] != 'reaction': # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    print(f"Final result: {boc_deprotection_at_late_stage}")
    return boc_deprotection_at_late_stage, findings_json
