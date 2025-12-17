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


ALCOHOL_FGS = [
    "Primary alcohol",
    "Secondary alcohol",
    "Tertiary alcohol",
    "Phenol",
]

HYDROXYL_PROTECTION_REACTIONS = [
    "Alcohol protection with silyl ethers",
    "Diol acetalization",
]

HYDROXYL_DEPROTECTION_REACTIONS = [
    "Alcohol deprotection from silyl ethers",
    "Alcohol deprotection from silyl ethers (double)",
    "Alcohol deprotection from silyl ethers (diol)",
    "Acetal hydrolysis to diol",
    "Cleavage of methoxy ethers to alcohols",
    "Cleavage of alkoxy ethers to alcohols",
    "Ether cleavage to primary alcohol",
    "Hydroxyl benzyl deprotection",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects hydroxyl protection or deprotection steps by identifying reactions from predefined lists.
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

    hydroxyl_protections = 0

    def dfs_traverse(node, depth, max_depth):
        nonlocal hydroxyl_protections, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                # Check for hydroxyl protection reactions
                fg_found_in_reactants = False
                for fg in ALCOHOL_FGS:
                    for r in reactants:
                        if checker.check_fg(fg, r):
                            if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(fg)
                            fg_found_in_reactants = True

                if fg_found_in_reactants:
                    for rxn in HYDROXYL_PROTECTION_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            hydroxyl_protections += 1
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            break # Found a protection reaction, no need to check others for this node

                # Check for hydroxyl deprotection reactions
                for rxn in HYDROXYL_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        hydroxyl_protections += 1
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break # Found a deprotection reaction, no need to check others for this node

        for child in node.get("children", []):
            # Propagating context. Assuming depth increases with recursion.
            # The actual values are not used in this specific function's logic.
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth, max_depth)

    # Initial call with dummy/placeholder context values as they are not provided to main.
    dfs_traverse(route, 0, 0)
    
    result = hydroxyl_protections >= 1
    
    if result:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "hydroxyl_protection_or_deprotection",
                "operator": ">=",
                "value": 1
            }
        })

    return result, findings_json