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


# Refactoring for Enumeration: Isolate lists of chemical entities
LATE_STAGE_PROTECTION_REACTIONS = [
    "Alcohol protection with silyl ethers",
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
    "Protection of carboxylic acid",
]

PROTECTING_GROUPS_OF_INTEREST = [
    "TMS ether protective group",
    "Silyl protective group",
    "Acetal/Ketal",
    "Boc",
    "Carbamic ester",  # Includes Cbz, Fmoc, etc.
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route involves a late-stage protection reaction (in the final 3 steps). It operates by checking for either specific, known protection reaction types or the formation of common protecting groups. The reaction types checked are defined in the `LATE_STAGE_PROTECTION_REACTIONS` list. The protecting groups checked are defined in the `PROTECTING_GROUPS_OF_INTEREST` list, and their formation is confirmed by their presence in the product and absence in all reactants.
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

    late_protection_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal late_protection_detected, findings_json

        if late_protection_detected:
            return

        if node["type"] == "reaction" and depth <= 2:  # Final 3 steps (depth 0, 1, 2)
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Check if this is a known protection reaction by name
                    for reaction_type in LATE_STAGE_PROTECTION_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            late_protection_detected = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": [
                                        "Alcohol protection with silyl ethers",
                                        "Boc amine protection",
                                        "Boc amine protection explicit",
                                        "Boc amine protection with Boc anhydride",
                                        "Boc amine protection (ethyl Boc)",
                                        "Boc amine protection of secondary amine",
                                        "Boc amine protection of primary amine",
                                        "Protection of carboxylic acid",
                                        "protecting_group_formation"
                                    ],
                                    "position": "final_3_steps"
                                }
                            })
                            return

                    # Check for the formation of common protecting groups
                    for pg in PROTECTING_GROUPS_OF_INTEREST:
                        # Check if protection group is in product and absent from ALL reactants
                        if checker.check_fg(pg, product_smiles) and not any(
                            checker.check_fg(pg, r) for r in reactants_smiles
                        ):
                            late_protection_detected = True
                            findings_json["atomic_checks"]["functional_groups"].append(pg)
                            findings_json["atomic_checks"]["named_reactions"].append("protecting_group_formation")
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": [
                                        "Alcohol protection with silyl ethers",
                                        "Boc amine protection",
                                        "Boc amine protection explicit",
                                        "Boc amine protection with Boc anhydride",
                                        "Boc amine protection (ethyl Boc)",
                                        "Boc amine protection of secondary amine",
                                        "Boc amine protection of primary amine",
                                        "Protection of carboxylic acid",
                                        "protecting_group_formation"
                                    ],
                                    "position": "final_3_steps"
                                }
                            })
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "protecting_group_formation",
                                        "PROTECTING_GROUPS_OF_INTEREST"
                                    ]
                                }
                            })
                            return
            except Exception:
                # Silently ignore errors in reaction analysis to prevent crashing
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return late_protection_detected, findings_json
