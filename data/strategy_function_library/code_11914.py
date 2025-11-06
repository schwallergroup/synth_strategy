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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthesis strategy that preserves the isoxazole heterocyclic core
    throughout the synthesis.

    In a retrosynthetic analysis, this checks if the isoxazole core in the final product
    is derived from a starting material that already contains this core, rather than
    being constructed during the synthesis.
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

    # Track if we found the heterocyclic core in the final product
    final_product_has_isoxazole = False
    result = False

    # Check if the final product has the required heterocycle
    if route["type"] == "mol":
        if checker.check_ring("isoxazole", route["smiles"]):
            final_product_has_isoxazole = True
            findings_json["atomic_checks"]["ring_systems"].append("isoxazole")
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "isoxazole",
                    "position": "final_product"
                }
            })

    # If final product doesn't have the core, return False immediately
    if not final_product_has_isoxazole:
        return False, findings_json

    # Track preservation through reactions
    core_preserved = True

    def check_preservation_across_reaction(reaction_node):
        nonlocal core_preserved, findings_json

        # Extract product and reactants from reaction SMILES
        try:
            rsmi = reaction_node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if any reactant has isoxazole
            reactants_with_isoxazole = [
                checker.check_ring("isoxazole", reactant) for reactant in reactants_smiles
            ]

            # Check if product has isoxazole
            product_has_isoxazole = checker.check_ring("isoxazole", product_smiles)
            if product_has_isoxazole and "isoxazole" not in findings_json["atomic_checks"]["ring_systems"]:
                findings_json["atomic_checks"]["ring_systems"].append("isoxazole")

            # Check preservation logic in retrosynthetic direction
            # If product has a core, at least one reactant should have it too
            if product_has_isoxazole and not any(reactants_with_isoxazole):
                core_preserved = False
                # If core is not preserved, it implies ring formation occurred
                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                # This also means the negation constraint is violated, so we don't add it.

        except (KeyError, Exception):
            # Any error during analysis means we can't confirm preservation
            core_preserved = False
            # If an error occurs and we can't confirm preservation, we assume ring formation might have happened
            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

    def dfs_traverse(node, depth=0):
        nonlocal core_preserved, findings_json
        # Process reaction nodes to check preservation
        if node["type"] == "reaction":
            check_preservation_across_reaction(node)

        # Determine the depth for the next level
        next_depth = depth
        if node["type"] != "reaction": # Only increase depth when going from chemical to reaction
            next_depth = depth + 1

        # Traverse children (retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Final check for the overall result
    result = core_preserved

    # Add negation constraint if core was preserved throughout
    if core_preserved:
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "ring_formation",
                "scope": "isoxazole"
            }
        })

    # Return True only if the core is found in final product and preserved throughout
    return result, findings_json
