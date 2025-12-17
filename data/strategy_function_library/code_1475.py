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
    This function detects a synthetic strategy where a nitrile group is preserved throughout the synthesis.
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

    nitrile_present_in_final = False
    all_steps_preserve_nitrile = True

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_present_in_final, all_steps_preserve_nitrile, findings_json

        if node["type"] == "mol":
            # Check if this is the final product (root node)
            if depth == 0:
                if checker.check_fg("Nitrile", node["smiles"]):
                    nitrile_present_in_final = True
                    if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

        elif node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has nitrile
                product_has_nitrile = checker.check_fg("Nitrile", product)
                if product_has_nitrile and "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

                # Check if any reactant has nitrile
                reactant_has_nitrile = any(checker.check_fg("Nitrile", r) for r in reactants)
                if reactant_has_nitrile and "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

                # If a reactant has nitrile but product doesn't, nitrile wasn't preserved
                if reactant_has_nitrile and not product_has_nitrile:
                    all_steps_preserve_nitrile = False

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Determine the new depth for recursive calls
        new_depth = depth
        if node["type"] != "reaction": # This means it's a 'mol' node, or any other type that should increase depth
            new_depth = depth + 1

        # Traverse children (retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    result = nitrile_present_in_final and all_steps_preserve_nitrile

    # Add structural constraints to findings_json based on the final result
    if nitrile_present_in_final:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Nitrile",
                "position": "last_stage"
            }
        })
    if all_steps_preserve_nitrile:
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "functional_group_consumption",
                "name": "Nitrile"
            }
        })

    return result, findings_json
