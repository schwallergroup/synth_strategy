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
    This function detects a synthetic strategy involving early-stage sulfonamide formation
    followed by ring transformations (opening and closing).
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

    # Track if we found the key features
    found_sulfonamide_formation = False
    found_ring_transformations = False

    def dfs_traverse(node, depth=0):
        nonlocal found_sulfonamide_formation, found_ring_transformations, findings_json

        # Update node depth
        node["depth"] = depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for sulfonamide formation (early stage, depth >= 2)
            if depth >= 2:
                # Check for sulfonamide formation reactions
                if checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                ):
                    print(f"Found early-stage sulfonamide formation at depth {depth}")
                    found_sulfonamide_formation = True
                    findings_json["atomic_checks"]["named_reactions"].append("Sulfonamide synthesis (Schotten-Baumann) primary amine")
                    if {"type": "positional", "details": {"target": "sulfonamide_formation", "position": "early_stage", "condition": "depth >= 2"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "sulfonamide_formation", "position": "early_stage", "condition": "depth >= 2"}})
                elif checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                ):
                    print(f"Found early-stage sulfonamide formation at depth {depth}")
                    found_sulfonamide_formation = True
                    findings_json["atomic_checks"]["named_reactions"].append("Sulfonamide synthesis (Schotten-Baumann) secondary amine")
                    if {"type": "positional", "details": {"target": "sulfonamide_formation", "position": "early_stage", "condition": "depth >= 2"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "sulfonamide_formation", "position": "early_stage", "condition": "depth >= 2"}})

                # Alternative check: product has sulfonamide but reactants don't
                elif checker.check_fg("Sulfonamide", product_smiles):
                    findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")
                    # Check if reactants don't have sulfonamide
                    if not any(checker.check_fg("Sulfonamide", r) for r in reactants_smiles if r):
                        print(f"Found early-stage sulfonamide formation at depth {depth}")
                        found_sulfonamide_formation = True
                        if {"type": "positional", "details": {"target": "sulfonamide_formation", "position": "early_stage", "condition": "depth >= 2"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "sulfonamide_formation", "position": "early_stage", "condition": "depth >= 2"}})

            # Check for ring transformations (late stage, depth 0 or 1)
            if depth <= 1:
                # Check if this is explicitly marked as a ring transformation
                if node["metadata"].get("RingBreaker", False):
                    print(f"Found ring transformation at depth {depth} (marked as RingBreaker)")
                    found_ring_transformations = True
                    findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")
                    if {"type": "positional", "details": {"target": "ring_destruction", "position": "late_stage", "condition": "depth <= 1"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ring_destruction", "position": "late_stage", "condition": "depth <= 1"}})

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, next is reaction
                new_depth = depth + 1
            # If current node is reaction, next is chemical, so depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Add co-occurrence constraint if both conditions are met
    result = found_sulfonamide_formation and found_ring_transformations
    if result:
        if {"type": "co-occurrence", "details": {"targets": ["sulfonamide_formation", "ring_destruction"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["sulfonamide_formation", "ring_destruction"]}})

    # Return True if both key features were found
    return result, findings_json
