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
    This function detects if a commercially available (in-stock) starting material in the synthetic route contains a trifluoromethyl group.
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

    has_cf3_building_block = False
    cf3_found_in_starting_material = False # Flag to track the specific structural constraint

    def dfs_traverse(node, depth=0):
        nonlocal has_cf3_building_block, findings_json, cf3_found_in_starting_material

        # Check molecules for CF3 groups
        if node["type"] == "mol":
            # Check if this is a starting material
            if node.get("in_stock", False):
                try:
                    mol_smiles = node["smiles"]
                    if checker.check_fg("Trifluoro group", mol_smiles):
                        has_cf3_building_block = True
                        findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
                        # If a trifluoro group is found in an in-stock material, set the flag for the structural constraint
                        cf3_found_in_starting_material = True
                except Exception as e:
                    pass

        # Traverse children
        # An early exit could be placed here, but is not added to adhere to modification constraints
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # After traversal, check if the structural constraint was met
    if cf3_found_in_starting_material:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Trifluoro group",
                    "starting_material"
                ],
                "scope": "single_molecule",
                "description": "A starting material (in-stock molecule) must contain a Trifluoro group."
            }
        })

    return has_cf3_building_block, findings_json
