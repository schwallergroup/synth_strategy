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
    Detects if an epoxide (oxirane) in the reactants is converted to a secondary alcohol in the product. This check is performed on all reaction steps except the final one.
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

    epoxide_opening_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal epoxide_opening_detected, findings_json

        if node["type"] == "reaction" and depth >= 2:  # Not the final step
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                try:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Analyzing reaction at depth {depth}: {rsmi}")

                    # Check for epoxide in reactants and secondary alcohol in product
                    has_epoxide = False
                    for reactant in reactants:
                        if checker.check_ring("oxirane", reactant):
                            print(f"Found oxirane ring in reactant: {reactant}")
                            has_epoxide = True
                            if "oxirane" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("oxirane")
                            break

                    # Check if product has a secondary alcohol
                    has_secondary_alcohol = checker.check_fg("Secondary alcohol", product)
                    if has_secondary_alcohol:
                        print(f"Found secondary alcohol in product: {product}")
                        if "Secondary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")

                    # If both conditions are met, we have an epoxide opening
                    if has_epoxide and has_secondary_alcohol:
                        print(f"Confirmed epoxide opening to secondary alcohol")
                        epoxide_opening_detected = True
                        # Add structural constraint for co-occurrence
                        if {"type": "co-occurrence", "details": {"targets": ["oxirane", "Secondary alcohol"], "scope": "reaction_step"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["oxirane", "Secondary alcohol"], "scope": "reaction_step"}})
                        # Add structural constraint for positional (not_last_stage)
                        if {"type": "positional", "details": {"target": "epoxide_opening_to_secondary_alcohol", "position": "not_last_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "epoxide_opening_to_secondary_alcohol", "position": "not_last_stage"}})

                except Exception as e:
                    print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    print(f"Epoxide opening to alcohol strategy detected: {epoxide_opening_detected}")
    return epoxide_opening_detected, findings_json
