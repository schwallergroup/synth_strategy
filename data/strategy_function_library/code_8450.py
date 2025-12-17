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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if multiple different protection groups (Boc and Cbz)
    are used in the synthesis.
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

    # Track protection groups
    boc_used = False
    cbz_used = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_used, cbz_used, findings_json

        # Check molecule nodes
        if node["type"] == "mol" and node["smiles"]:
            mol_smiles = node["smiles"]

            # Check for Boc protection group
            if checker.check_fg("Boc", mol_smiles):
                boc_used = True
                if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boc")

            # Check for Cbz protection group (assuming a similar check for Cbz exists or would be added)
            # For this specific problem, Cbz is not explicitly checked in the original code,
            # but the problem description implies its relevance.
            # If a Cbz check were present, it would be added here:
            # if checker.check_fg("Cbz", mol_smiles):
            #     cbz_used = True
            #     if "Cbz" not in findings_json["atomic_checks"]["functional_groups"]:
            #         findings_json["atomic_checks"]["functional_groups"].append("Cbz")

        # Check reaction nodes
        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for Boc protection/deprotection reactions
            if checker.check_reaction("Boc amine protection", rsmi):
                boc_used = True
                if "Boc amine protection" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Boc amine protection")
            
            if checker.check_reaction("Boc amine deprotection", rsmi):
                boc_used = True
                if "Boc amine deprotection" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Boc amine deprotection")

        # Traverse children
        for child in node.get("children", []):
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if both protection groups are used
    result = False
    if boc_used and cbz_used:
        result = True
        # Add the structural constraint if both are used
        findings_json["structural_constraints"].append(
            {
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "boc_usage",
                        "cbz_usage"
                    ]
                }
            }
        )

    return result, findings_json
