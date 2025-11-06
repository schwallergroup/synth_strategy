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
    This function detects if an aromatic scaffold (like chloronaphthalene)
    is preserved throughout the synthesis.
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

    # Track reactions where chloronaphthalene scaffold is preserved
    preserved_reactions = []
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal preserved_reactions, findings_json
        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants = parts[0].split(".")
            product = parts[-1]

            # Check if chloronaphthalene scaffold exists in product
            product_has_scaffold = (
                checker.check_ring("naphthalene", product)
                and "Cl" in product
            )
            if checker.check_ring("naphthalene", product):
                if "naphthalene" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("naphthalene")
            if "Cl" in product:
                if "chloro" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("chloro")

            # Check if chloronaphthalene scaffold exists in any reactant
            reactant_has_scaffold = False
            for reactant in reactants:
                if (
                    checker.check_ring("naphthalene", reactant)
                    and "Cl" in reactant
                ):
                    reactant_has_scaffold = True
                    if "naphthalene" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("naphthalene")
                    if "chloro" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("chloro")
                    break

            # If scaffold exists in both reactant and product, it's preserved in this reaction
            if reactant_has_scaffold and product_has_scaffold:
                preserved_reactions.append((depth, rsmi))
                if "scaffold_preservation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("scaffold_preservation")
                print(f"Scaffold preserved at depth {depth}")

        # Determine the new depth for the recursive call
        new_depth = depth
        if node["type"] != "reaction":  # This means it's a chemical node
            new_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth (in retrosynthetic direction)
    preserved_reactions.sort(key=lambda x: x[0])

    # Check for 3 consecutive reactions with preserved scaffold
    if len(preserved_reactions) >= 3:
        # Find sequences of consecutive depths
        consecutive_count = 1
        max_consecutive = 1

        for i in range(1, len(preserved_reactions)):
            current_depth = preserved_reactions[i][0]
            prev_depth = preserved_reactions[i - 1][0]

            # Check if depths are consecutive (may have gaps of 2 in retrosynthetic direction)
            if current_depth - prev_depth <= 2:
                consecutive_count += 1
                max_consecutive = max(max_consecutive, consecutive_count)
            else:
                consecutive_count = 1

        if max_consecutive >= 3:
            print(
                f"Chloronaphthalene scaffold preserved across {max_consecutive} consecutive reactions"
            )
            result = True
            # Add structural constraint finding
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "scaffold_preservation",
                    "operator": ">=",
                    "value": 3,
                    "condition": "consecutive_steps",
                    "consecutive_details": {
                        "description": "Counts the maximum number of scaffold preservation events that occur in a sequence where the retrosynthetic depth between adjacent events is no more than 2.",
                        "max_depth_difference": 2
                    }
                }
            })

    return result, findings_json
