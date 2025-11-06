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


def main(route, max_depth) -> Tuple[bool, Dict]:
    """
    Detects a strategy where a ketone, present in an intermediate, is used as a functional handle
    for a late-stage reduction to a secondary alcohol.
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

    # Track if we found ketone being used as a handle
    found_ketone_to_alcohol = False
    ketone_present_in_middle = False

    def dfs_traverse(node, max_depth, depth=0):
        nonlocal found_ketone_to_alcohol, ketone_present_in_middle, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            product_smiles = rsmi.split(">")[-1]

            # Check for ketone to alcohol transformation in late stage (final or penultimate step)
            if depth <= 1:
                # Check if this is a ketone reduction reaction
                if checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi):
                    print(f"Found ketone reduction to alcohol at depth {depth}")
                    found_ketone_to_alcohol = True
                    if "Reduction of ketone to secondary alcohol" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Reduction of ketone to secondary alcohol")
                    if {"type": "positional", "details": {"target": "Reduction of ketone to secondary alcohol", "position": "last_or_penultimate_stage"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Reduction of ketone to secondary alcohol", "position": "last_or_penultimate_stage"}})

            # Check for ketone presence in middle stages
            if 2 <= depth <= 3:
                if checker.check_fg("Ketone", product_smiles):
                    print(f"Found ketone present in middle stage at depth {depth}")
                    ketone_present_in_middle = True
                    if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                    if {"type": "positional", "details": {"target": "Ketone", "position": "middle_stage"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Ketone", "position": "middle_stage"}})

        elif node["type"] == "mol":
            # Also check for ketones in molecule nodes in middle stages
            if 2 <= depth <= 3:
                if checker.check_fg("Ketone", node["smiles"]):
                    print(f"Found ketone in molecule at depth {depth}")
                    ketone_present_in_middle = True
                    if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                    if {"type": "positional", "details": {"target": "Ketone", "position": "middle_stage"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Ketone", "position": "middle_stage"}})

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, max_depth, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, max_depth, depth + 1)

    # Start traversal
    dfs_traverse(route, max_depth, 0)

    print(f"Ketone to alcohol found: {found_ketone_to_alcohol}")
    print(f"Ketone in middle stages: {ketone_present_in_middle}")

    # Return True if ketone was used as a functional handle
    result = found_ketone_to_alcohol and ketone_present_in_middle

    if result:
        if {"type": "co-occurrence", "details": {"targets": ["Reduction of ketone to secondary alcohol", "Ketone"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Reduction of ketone to secondary alcohol", "Ketone"]}})

    return result, findings_json
