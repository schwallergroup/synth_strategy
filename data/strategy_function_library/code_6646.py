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
    This function detects a synthetic strategy involving late-stage cyanation
    (introduction of a nitrile group in the final step).
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

    final_reaction_has_cyanation = False
    earlier_reactions_have_cyanation = False

    def dfs_traverse(node, depth=0):
        nonlocal final_reaction_has_cyanation
        nonlocal earlier_reactions_have_cyanation
        nonlocal findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has a nitrile group
                product_has_nitrile = checker.check_fg("Nitrile", product)
                if product_has_nitrile:
                    if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

                # A reaction is a cyanation if a nitrile is formed.
                is_cyanation = False
                if product_has_nitrile:
                    reactant_has_nitrile = False
                    for r in reactants:
                        if checker.check_fg("Nitrile", r):
                            reactant_has_nitrile = True
                            break
                    if not reactant_has_nitrile:
                        is_cyanation = True
                        if "cyanation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("cyanation")

                if depth == 1 and is_cyanation:
                    final_reaction_has_cyanation = True
                elif depth > 1 and is_cyanation:
                    earlier_reactions_have_cyanation = True
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)

    # True only if cyanation occurs in final step but not in earlier steps
    result = final_reaction_has_cyanation and not earlier_reactions_have_cyanation

    if final_reaction_has_cyanation:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "cyanation",
                "position": "last_stage"
            }
        })
    if earlier_reactions_have_cyanation:
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "cyanation",
                "position": "not_last_stage"
            }
        })

    return result, findings_json