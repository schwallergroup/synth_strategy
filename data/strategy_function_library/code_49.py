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
    This function detects Boc deprotection steps by checking if the reaction matches any of the names in a predefined list of deprotection reactions: 'Boc amine deprotection', 'Boc amine deprotection of guanidine', 'Boc amine deprotection to NH-NH2', and 'Tert-butyl deprotection of amine'.
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

    boc_deprotection = False
    # This assumes max_depth can be determined from the route object, a common pattern.
    max_depth = route.get("max_depth", 0)

    def dfs_traverse(node, reaction, depth, max_depth):
        nonlocal boc_deprotection, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for Boc deprotection reactions using the checker function and the refactored list
            for rxn_name in BOC_DEPROTECTION_REACTIONS:
                if checker.check_reaction(rxn_name, rsmi):
                    print(f"Boc deprotection detected in reaction: {rsmi}")
                    boc_deprotection = True
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                    # If multiple reactions could match, we only need to record one for this specific check
                    # If we wanted all, we'd remove the break and ensure unique entries later.
                    break

        # Determine the new depth based on the current node's type
        new_depth = depth
        if node["type"] != "reaction": # This means it's a chemical node
            new_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            # In a real system, the reaction object for the child would be passed here.
            dfs_traverse(child, None, new_depth, max_depth)

    # Start traversal from root
    dfs_traverse(route, None, 1, max_depth)
    return boc_deprotection, findings_json
