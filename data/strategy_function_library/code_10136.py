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
    Detects synthesis routes where a nitrile group is present throughout the synthesis
    and not modified.
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

    result = True

    # Check if target molecule has nitrile
    if route["type"] == "mol":
        if not checker.check_fg("Nitrile", route["smiles"]):
            print("Target molecule does not have a nitrile group")
            result = False
        else:
            findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "Nitrile",
                    "position": "final_product"
                }
            })

    nitrile_preserved = True

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_preserved, result, findings_json

        if not nitrile_preserved:  # Early termination if already False
            return

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitrile in product
                product_has_nitrile = checker.check_fg("Nitrile", product)
                if product_has_nitrile and "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

                # Check if at least one reactant has nitrile
                reactant_has_nitrile = any(checker.check_fg("Nitrile", r) for r in reactants)
                if reactant_has_nitrile and "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

                # If product has nitrile but no reactant has it, nitrile was created (not preserved)
                if product_has_nitrile and not reactant_has_nitrile:
                    print(f"Nitrile was created in reaction, not preserved: {rsmi}")
                    nitrile_preserved = False
                    result = False
                    if "functional_group_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("functional_group_formation")
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "functional_group_formation",
                            "scope": "Nitrile"
                        }
                    })

                # If reactant has nitrile but product doesn't, nitrile was modified
                if reactant_has_nitrile and not product_has_nitrile:
                    print(f"Nitrile was modified in reaction: {rsmi}")
                    nitrile_preserved = False
                    result = False
                    if "functional_group_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("functional_group_destruction")
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "functional_group_destruction",
                            "scope": "Nitrile"
                        }
                    })

            except Exception as e:
                print(f"Error processing reaction: {e}")
                nitrile_preserved = False
                result = False

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., chemical)
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # The final result is determined by nitrile_preserved after traversal
    if not nitrile_preserved:
        result = False

    print(f"Nitrile preserved throughout synthesis: {result}")
    return result, findings_json