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
    Detects SEM protection-deprotection strategy.
    Looks for SEM group in intermediates and its removal in a late stage.
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

    sem_protected_intermediate = False
    sem_deprotection_at_late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal sem_protected_intermediate, sem_deprotection_at_late_stage, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for SEM protection in reactants
                for reactant in reactants:
                    if checker.check_fg("SEM", reactant):
                        sem_protected_intermediate = True
                        if "SEM" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("SEM")
                        # print(f"Found SEM-protected intermediate at depth {depth}")

                # Check for SEM deprotection by its presence in reactants and absence in product
                reactant_has_sem = any(
                    checker.check_fg("SEM", r) for r in reactants
                )
                product_has_sem = checker.check_fg("SEM", product)

                if (reactant_has_sem and not product_has_sem):
                    if "SEM" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("SEM")
                    # Assuming "SEM_deprotection" is a named reaction if this condition is met
                    if "SEM_deprotection" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("SEM_deprotection")

                    if depth <= 1:
                        sem_deprotection_at_late_stage = True
                        # print(f"Found late-stage SEM deprotection at depth {depth}")

        # For molecule nodes, check if they contain SEM groups
        elif node["type"] == "mol" and not node.get("in_stock", False):
            if checker.check_fg("SEM", node["smiles"]):
                sem_protected_intermediate = True
                if "SEM" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("SEM")
                # print(f"Found SEM-protected molecule at depth {depth}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = sem_protected_intermediate and sem_deprotection_at_late_stage

    if result:
        # Add structural constraints if both conditions are met
        findings_json["structural_constraints"].append(
            {
                "type": "positional",
                "details": {
                    "target": "SEM_deprotection",
                    "position": "last_two_stages"
                }
            }
        )
        findings_json["structural_constraints"].append(
            {
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "SEM",
                        "SEM_deprotection"
                    ]
                }
            }
        )

    return result, findings_json
