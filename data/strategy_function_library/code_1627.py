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
    This function detects if the synthesis involves late-stage isoxazole formation
    as the final step in a multi-ring system synthesis.
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
    depth_of_formation = None
    result = False

    def count_rings(mol_smiles):
        """Count the number of rings in a molecule"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return 0
        return mol.GetRingInfo().NumRings()

    def dfs_traverse(node, depth=0):
        nonlocal isoxazole_formed, depth_of_formation, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains isoxazole
            has_isoxazole_in_product = checker.check_ring("isoxazole", product)

            if has_isoxazole_in_product:
                if "isoxazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("isoxazole")

                # Check if any reactant contains isoxazole
                reactants_have_isoxazole = any(
                    checker.check_ring("isoxazole", reactant)
                    for reactant in reactants
                )

                if not reactants_have_isoxazole:
                    # Isoxazole formation detected, now check for multi-ring system
                    product_ring_count = count_rings(product)

                    if product_ring_count > 1:
                        isoxazole_formed = True
                        depth_of_formation = depth
                        # Record the count constraint if met
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "rings_in_product_of_isoxazole_formation",
                                "operator": ">",
                                "value": 1
                            }
                        })

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if isoxazole formation is late-stage (depth = 1)
    if isoxazole_formed and depth_of_formation == 1:
        result = True
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "isoxazole_formation",
                "position": "last_stage"
            }
        })

    # Also check for depth 0 as it might be the final step in some cases
    if isoxazole_formed and depth_of_formation == 0:
        result = True
        # Only add if not already added by depth 1 check (to avoid duplicates if both conditions are met)
        # In this specific case, depth 0 and 1 are mutually exclusive for a single formation event,
        # but good practice to check for duplicates if the constraint could be met by multiple paths.
        if {"type": "positional", "details": {"target": "isoxazole_formation", "position": "last_stage"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "isoxazole_formation",
                    "position": "last_stage"
                }
            })

    return result, findings_json
