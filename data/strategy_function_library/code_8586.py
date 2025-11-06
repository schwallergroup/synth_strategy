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
    This function detects late-stage isoxazole formation strategy.
    It looks for the creation of an isoxazole ring in the final or penultimate synthetic step.
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

    isoxazole_formed = False

    def dfs_traverse(node, depth=0):
        nonlocal isoxazole_formed, findings_json

        # Check final step (depth 0) and penultimate step (depth 1)
        if node["type"] == "reaction" and depth <= 1:
            # Extract reactants and product
            try:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    return

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if isoxazole is in product
                if checker.check_ring("isoxazole", product_smiles):
                    findings_json["atomic_checks"]["ring_systems"].append("isoxazole")

                    # Check if isoxazole is not in any reactant
                    isoxazole_in_reactants = False
                    for reactant in reactants_smiles:
                        if checker.check_ring("isoxazole", reactant):
                            isoxazole_in_reactants = True
                            break

                    if not isoxazole_in_reactants:
                        isoxazole_formed = True
                        # Record the structural constraint if the condition is met
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "isoxazole_ring_formation",
                                "position": "late_stage"
                            }
                        })
                        # Assuming 'ring_formation' is a general reaction type that applies here
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

            except Exception as e:
                pass

        # Continue traversal
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)
    print(f"Final result: isoxazole_formed = {isoxazole_formed}")

    return isoxazole_formed, findings_json
