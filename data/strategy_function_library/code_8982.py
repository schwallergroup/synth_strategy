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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthesis strategy where both trifluoromethyl group and
    oxazole heterocycle are preserved throughout the synthesis.
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

    result = True

    # Check if target molecule contains the groups
    target_has_cf3 = checker.check_fg("Trifluoro group", route["smiles"])
    if target_has_cf3:
        findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
    target_has_oxazole = checker.check_ring("oxazole", route["smiles"])
    if target_has_oxazole:
        findings_json["atomic_checks"]["ring_systems"].append("oxazole")

    # If target doesn't have both groups, return False
    if not (target_has_cf3 and target_has_oxazole):
        print(f"Target molecule doesn't contain both trifluoromethyl and oxazole groups")
        result = False
    else:
        # Add structural constraint for co-occurrence in final product
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Trifluoro group",
                    "oxazole"
                ],
                "scope": "final_product"
            }
        })

    # Track if each reaction preserves these groups
    all_steps_preserve_cf3 = True
    all_steps_preserve_oxazole = True

    def dfs_traverse(node, depth=0):
        nonlocal all_steps_preserve_cf3, all_steps_preserve_oxazole, result, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product contains the groups
            product_has_cf3 = checker.check_fg("Trifluoro group", product_smiles)
            product_has_oxazole = checker.check_ring("oxazole", product_smiles)

            # Check if any reactant contains the groups
            reactants_have_cf3 = any(
                checker.check_fg("Trifluoro group", r) for r in reactants_smiles
            )
            reactants_have_oxazole = any(checker.check_ring("oxazole", r) for r in reactants_smiles)

            # In retrosynthesis, we check if a group in the product is preserved in at least one reactant
            if product_has_cf3 and not reactants_have_cf3:
                all_steps_preserve_cf3 = False
                print(f"Trifluoromethyl group not preserved at depth {depth}")

            if product_has_oxazole and not reactants_have_oxazole:
                all_steps_preserve_oxazole = False
                print(f"Oxazole heterocycle not preserved at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Update result based on preservation checks
    if not all_steps_preserve_cf3:
        result = False
    else:
        # Add structural constraint for negation of formation of Trifluoro group
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "formation of Trifluoro group",
                "scope": "all_steps"
            }
        })

    if not all_steps_preserve_oxazole:
        result = False
    else:
        # Add structural constraint for negation of formation of oxazole
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "formation of oxazole",
                "scope": "all_steps"
            }
        })

    # Return True if both groups are preserved throughout
    return result, findings_json
