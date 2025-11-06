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
    This function detects preservation of CF3 groups throughout the synthesis.

    It checks if CF3 groups are present in the target molecule and preserved
    throughout all synthesis steps.
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

    # Track if CF3 is preserved throughout synthesis
    cf3_preserved = True

    # Track if target molecule has CF3
    target_has_cf3 = False

    # Track molecules with CF3 by their SMILES
    molecules_with_cf3 = set()

    def dfs_traverse(node, depth=0):
        nonlocal cf3_preserved, target_has_cf3, molecules_with_cf3, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check for CF3 group using the provided checker function
            has_cf3 = checker.check_fg("Trifluoro group", mol_smiles)

            if has_cf3:
                molecules_with_cf3.add(mol_smiles)
                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

            # If this is the target molecule (root node)
            if depth == 0:
                target_has_cf3 = has_cf3
                print(f"Target molecule: {mol_smiles}")
                print(f"Target has CF3 group: {has_cf3}")
                if has_cf3:
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Trifluoro group",
                            "position": "last_stage"
                        }
                    })

                # If target doesn't have CF3, no need to check preservation
                if not has_cf3:
                    cf3_preserved = False
                    return

            # For starting materials, check if they contribute CF3
            elif node.get("in_stock", False) and has_cf3:
                print(f"Starting material with CF3 found: {mol_smiles}")

        # Traverse children
        children_with_cf3 = False
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)
            if child["type"] == "mol" and child["smiles"] in molecules_with_cf3:
                children_with_cf3 = True

        # After traversing all children of a reaction node, check if CF3 was preserved
        if node["type"] == "reaction" and target_has_cf3:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                product = rsmi.split(">")[-1]

                # If any reactant (child) has CF3 but the product does not, preservation has failed.
                if children_with_cf3 and not checker.check_fg("Trifluoro group", product):
                    cf3_preserved = False
                    print(f"CF3 group lost in reaction (retrosynthetic direction): {rsmi}")
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "description": "A reaction step must not result in the loss of a Trifluoro group. A loss is defined as a reaction where at least one reactant contains a Trifluoro group, but the product does not."
                        }
                    })
            except Exception as e:
                print(f"Error in post-children analysis: {e}")

    # Start traversal from root
    dfs_traverse(route)

    # Final determination
    result = False
    if target_has_cf3:
        # Check if CF3 was preserved
        if cf3_preserved:
            print("CF3 group preserved throughout the synthesis")
            result = True
        else:
            print("CF3 group not preserved in all synthesis steps")
            result = False
    else:
        print("Target molecule does not contain a CF3 group")
        result = False

    return result, findings_json
