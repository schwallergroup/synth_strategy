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
    Detects if a purine scaffold is introduced in the late stage of the synthesis (low depth).
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

    purine_introduced = False
    purine_modified = False
    purine_depth = float("inf")
    result = False

    # First check if the target molecule contains purine
    if route["type"] == "mol" and route["smiles"]:
        target_mol_smiles = route["smiles"]
        has_purine = checker.check_ring("purine", target_mol_smiles)
        print(f"Target molecule has purine: {has_purine}")
        if has_purine:
            findings_json["atomic_checks"]["ring_systems"].append("purine")
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "purine"
                    ],
                    "scope": "final_product"
                }
            })
        if not has_purine:
            print("Target molecule does not contain purine, returning False")
            return False, findings_json

    def dfs_traverse(node, depth=0):
        nonlocal purine_introduced, purine_modified, purine_depth, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            try:
                # Check if product contains purine
                product_has_purine = checker.check_ring("purine", product)
                print(f"Product has purine: {product_has_purine}")

                if product_has_purine:
                    # Check if any reactant contains purine
                    reactants_with_purine = [
                        r for r in reactants if checker.check_ring("purine", r)
                    ]
                    reactants_have_purine = len(reactants_with_purine) > 0

                    for reactant in reactants_with_purine:
                        print(f"Reactant has purine: {reactant}")

                    print(f"Reactants have purine: {reactants_have_purine}")

                    # If product has purine but reactants don't, purine is introduced in this step
                    if not reactants_have_purine:
                        purine_introduced = True
                        purine_depth = min(purine_depth, depth)
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                        print(f"Purine introduced at depth {depth}: {rsmi}")
                    # If both product and reactants have purine, purine is being modified
                    elif reactants_have_purine:
                        purine_modified = True
                        purine_depth = min(purine_depth, depth)
                        findings_json["atomic_checks"]["named_reactions"].append("ring_modification")
                        print(f"Purine modified at depth {depth}: {rsmi}")
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Purine introduced or modified: {purine_introduced or purine_modified}, at depth: {purine_depth}"
    )

    # Return True if purine is introduced or modified in the late stage (depth <= 2)
    result = (purine_introduced or purine_modified) and purine_depth <= 2

    if result:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": [
                    "ring_formation",
                    "ring_modification"
                ],
                "position": "late_stage",
                "max_depth": 2
            }
        })

    return result, findings_json