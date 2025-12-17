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
    Detects if indole cores are preserved throughout the synthesis route.

    Returns True if indole cores are preserved in all non-starting material
    molecules throughout the synthesis route.
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

    # Track if we've found any indole preservation issues
    preservation_issue_found = False

    def dfs_traverse(node, depth=0):
        nonlocal preservation_issue_found, findings_json

        if preservation_issue_found:
            return  # Early termination if we already found an issue

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_indole = checker.check_ring("indole", mol_smiles)
            if has_indole and "indole" not in findings_json["atomic_checks"]["ring_systems"]:
                findings_json["atomic_checks"]["ring_systems"].append("indole")

            # If this is a starting material, no need to check further
            if node.get("in_stock", False):
                return

            # Continue traversal
            for child in node.get("children", []):
                if child["type"] == "reaction":
                    # Check if indole is preserved in this reaction
                    rsmi = child["metadata"].get("rsmi", "")
                    if rsmi:
                        product = rsmi.split(">")[-1]
                        reactants = rsmi.split(">")[0].split(".")

                        # In retrosynthesis: product is current molecule, reactants are precursors
                        if has_indole:
                            # If product has indole, at least one reactant should have indole
                            reactant_has_indole = any(
                                checker.check_ring("indole", reactant) for reactant in reactants
                            )

                            if not reactant_has_indole:
                                # Check if all reactants are starting materials
                                all_starting_materials = True
                                for grandchild in child.get("children", []):
                                    if grandchild["type"] == "mol" and not grandchild.get(
                                        "in_stock", False
                                    ):
                                        all_starting_materials = False
                                        break

                                if not all_starting_materials:
                                    preservation_issue_found = True
                                    # Record the negation constraint violation
                                    findings_json["structural_constraints"].append({
                                        "type": "negation",
                                        "details": {
                                            "description": "The strategy fails if an indole ring is formed in any reaction step that uses at least one non-starting-material precursor.",
                                            "event": "ring_formation",
                                            "scope": "indole",
                                            "positional_context": "not_first_stage"
                                        }
                                    })
                                    return

                    # Depth increases when going from chemical to reaction
                    dfs_traverse(child, depth + 1)
                else:
                    # Depth increases when going from chemical to chemical (shouldn't happen in this graph structure)
                    dfs_traverse(child, depth + 1)

        elif node["type"] == "reaction":
            # Continue traversal for all children
            for child in node.get("children", []):
                # Depth remains the same when going from reaction to chemical
                dfs_traverse(child, depth)

    result = False
    # Check if target molecule has indole
    if route["type"] == "mol" and checker.check_ring("indole", route["smiles"]):
        if "indole" not in findings_json["atomic_checks"]["ring_systems"]:
            findings_json["atomic_checks"]["ring_systems"].append("indole")
        # Record the positional constraint for the target molecule having indole
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "indole",
                "position": "last_stage"
            }
        })
        # Start traversal from target molecule
        dfs_traverse(route)
        # If no preservation issues found, strategy succeeds
        result = not preservation_issue_found
    else:
        # If target doesn't have indole, strategy doesn't apply
        result = False

    return result, findings_json
