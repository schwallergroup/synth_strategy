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
    Detects if a nitrile functional group is preserved throughout the entire synthesis.
    This means the target molecule has a nitrile group, and this group is neither created
    nor destroyed in any reaction step of the synthesis.
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

    # First check if the target molecule (root of the tree) has a nitrile group
    if not checker.check_fg("Nitrile", route["smiles"]):
        print(f"Target molecule does not contain a nitrile group: {route['smiles']}")
        result = False
    else:
        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Nitrile",
                "position": "last_stage"
            }
        })

    print(f"Target molecule contains nitrile: {route['smiles']}")

    # Track if nitrile is preserved in all reactions
    nitrile_preserved = True

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_preserved, findings_json, result

        if not nitrile_preserved:
            return  # Early termination

        # Process reaction nodes to check if nitrile is preserved
        if node["type"] == "reaction":
            try:
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

                # If product has nitrile but no reactant has nitrile, then nitrile was created
                if product_has_nitrile and not reactant_has_nitrile:
                    print(f"Nitrile was created in reaction: {rsmi}")
                    nitrile_preserved = False
                    result = False
                    if "Nitrile_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Nitrile_formation")
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "Nitrile_formation"
                        }
                    })

                # If any reactant has nitrile but product doesn't, then nitrile was destroyed
                if reactant_has_nitrile and not product_has_nitrile:
                    print(f"Nitrile was destroyed in reaction: {rsmi}")
                    nitrile_preserved = False
                    result = False
                    if "Nitrile_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Nitrile_destruction")
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "Nitrile_destruction"
                        }
                    })

            except Exception as e:
                print(f"Error analyzing reaction: {e}")
                nitrile_preserved = False
                result = False

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases for children
                new_depth = depth + 1
            # If current node is reaction, depth remains the same for children
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    if nitrile_preserved:
        print("Nitrile group is preserved throughout the synthesis")
    else:
        print("Nitrile group is not preserved throughout the synthesis")

    return result, findings_json
