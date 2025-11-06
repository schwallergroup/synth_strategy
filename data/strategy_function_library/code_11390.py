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
    This function detects a synthetic strategy involving the use of a
    trifluoromethyl-containing building block.
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

    def dfs_traverse(node, depth=0):
        nonlocal has_cf3_building_block, findings_json

        if node.get("type") == "mol":
            mol_smiles = node.get("smiles", "")
            if mol_smiles and checker.check_fg("Trifluoro group", mol_smiles):
                findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
                if node.get("in_stock", False) or depth >= 5:
                    has_cf3_building_block = True
                    # Add the structural constraint if the condition is met
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Trifluoro group",
                            "position": "on_starting_material"
                        }
                    })

        # Traverse children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            if node.get("type") == "reaction":
                new_depth = depth
            else:  # Assuming 'mol' or 'chemical' type, or any other type that should increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from root
    print(f"Starting analysis of synthetic route")
    dfs_traverse(route)
    print(f"Analysis complete. Found CF3 building block: {has_cf3_building_block}")
    return has_cf3_building_block, findings_json
