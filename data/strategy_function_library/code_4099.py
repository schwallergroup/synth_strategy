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
    This function detects a synthetic strategy involving functionalization of a thiophene heterocycle,
    particularly focusing on late-stage formylation of thiophene.
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

    # Track if we found thiophene functionalization
    found_thiophene_formylation = False
    max_depth = 0
    formylation_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal found_thiophene_formylation, max_depth, formylation_depth, findings_json

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    return

                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains thiophene
                if checker.check_ring("thiophene", product):
                    if "thiophene" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("thiophene")
                    # Check if product has aldehyde
                    if checker.check_fg("Aldehyde", product):
                        if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                        # Check if any reactant has thiophene
                        thiophene_in_reactants = any(
                            checker.check_ring("thiophene", r) for r in reactants
                        )

                        # Check if any reactant has aldehyde on thiophene
                        aldehyde_in_reactants = any(
                            checker.check_fg("Aldehyde", r) and checker.check_ring("thiophene", r)
                            for r in reactants
                        )

                        # If thiophene is in reactants but aldehyde is not on thiophene in reactants,
                        # then this is a formylation step
                        if thiophene_in_reactants and not aldehyde_in_reactants:
                            print(f"Found thiophene formylation at depth {depth}")
                            print(f"Reaction SMILES: {rsmi}")
                            found_thiophene_formylation = True
                            formylation_depth = min(formylation_depth, depth)
                            if "formylation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("formylation")
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if formylation is late-stage (in the first third of the synthesis)
    is_late_stage = formylation_depth < max_depth / 3 if max_depth > 0 else False

    print(f"Found thiophene formylation: {found_thiophene_formylation}")
    print(f"Formylation depth: {formylation_depth}, Max depth: {max_depth}")
    print(f"Is late-stage: {is_late_stage}")

    # Return True if thiophene was functionalized with formylation in a late stage
    result = found_thiophene_formylation and is_late_stage

    if result:
        # Add the structural constraint if the overall condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "formylation_on_thiophene",
                "position": "late_stage (top third of route)"
            }
        })

    return result, findings_json
