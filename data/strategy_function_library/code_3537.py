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
    This function detects a synthetic strategy where a fluorinated aryl group
    is maintained throughout the synthesis.
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

    # Track molecules with fluorinated aryl groups
    molecules_with_fluoro_aryl = []

    def dfs_traverse(node, depth=0):
        nonlocal molecules_with_fluoro_aryl, findings_json
        if node["type"] == "mol":
            # Skip empty SMILES
            if not node["smiles"]:
                return

            try:
                # Check if molecule contains fluorinated aryl group
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Use checker to detect 'aryl fluoride'
                    if checker.check_fg("aryl fluoride", mol):
                        molecules_with_fluoro_aryl.append((node["smiles"], depth))
                        if "aryl fluoride" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("aryl fluoride")
                        print(f"Found fluorinated aryl at depth {depth}: {node['smiles']}")
            except Exception as e:
                print(f"Error processing SMILES: {node['smiles']} - {str(e)}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is 'mol' (chemical), depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Sort molecules by depth (retrosynthetic order)
    molecules_with_fluoro_aryl.sort(key=lambda x: x[1])

    # Strategy is present if:
    # 1. At least 3 molecules in the route contain fluorinated aryl groups
    # 2. The final product (depth 0) contains a fluorinated aryl group
    strategy_present = False

    if len(molecules_with_fluoro_aryl) >= 3:
        # Check if final product (lowest depth) has fluorinated aryl
        final_product_has_fluoro = any(depth == 0 for _, depth in molecules_with_fluoro_aryl)

        # Check if at least 2 intermediates have fluorinated aryl
        intermediates_with_fluoro = sum(1 for _, depth in molecules_with_fluoro_aryl if depth > 0)

        strategy_present = final_product_has_fluoro and intermediates_with_fluoro >= 2

        if strategy_present:
            # Add structural constraints if the strategy is detected
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "aryl fluoride",
                    "operator": ">=",
                    "value": 3
                }
            })
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "aryl fluoride",
                    "position": "last_stage"
                }
            })

    print(f"Fluorinated aryl maintenance strategy detected: {strategy_present}")
    print(f"Total molecules with fluorinated aryl: {len(molecules_with_fluoro_aryl)}")

    return strategy_present, findings_json
