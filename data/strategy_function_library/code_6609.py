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
    This function detects a synthesis strategy where a chloro-aryl substituent is preserved
    throughout the synthesis.
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

    # Track molecule nodes with chloro-aryl
    mol_nodes_with_chloro_aryl = []
    mol_nodes_total = []

    def dfs_traverse(node, depth=0):
        nonlocal mol_nodes_with_chloro_aryl, mol_nodes_total, findings_json
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for chloro-aryl substituent using the checker function
            has_chloro_aryl = checker.check_fg("Aromatic chloride", mol_smiles)

            if has_chloro_aryl:
                mol_nodes_with_chloro_aryl.append((depth, mol_smiles))
                findings_json["atomic_checks"]["functional_groups"].append("Aromatic chloride")
                print(f"Found chloro-aryl substituent at depth {depth}")

            mol_nodes_total.append((depth, mol_smiles))

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # Only increase depth if not a reaction node
            next_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    result = False

    # No molecule nodes found
    if not mol_nodes_total:
        return result, findings_json

    # Check if we have at least 2 molecule nodes with chloro-aryl
    if len(mol_nodes_with_chloro_aryl) >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "molecules_with_Aromatic_chloride",
                "operator": ">=",
                "value": 2
            }
        })

    # Check if the final product (depth 0) has a chloro-aryl
    final_product_has_chloro_aryl = any(depth == 0 for depth, _ in mol_nodes_with_chloro_aryl)
    if final_product_has_chloro_aryl:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Aromatic chloride",
                "position": "last_stage"
            }
        })

    # Get the maximum depth to find starting materials
    max_depth = max(depth for depth, _ in mol_nodes_total) if mol_nodes_total else -1

    # Check if at least one starting material has a chloro-aryl
    starting_material_has_chloro_aryl = any(
        depth == max_depth for depth, _ in mol_nodes_with_chloro_aryl
    )
    if starting_material_has_chloro_aryl:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Aromatic chloride",
                "position": "first_stage"
            }
        })

    # Calculate the percentage of molecule nodes that have chloro-aryl
    preservation_ratio = len(mol_nodes_with_chloro_aryl) / len(mol_nodes_total) if mol_nodes_total else 0

    if preservation_ratio >= 0.5:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "ratio_of_molecules_with_Aromatic_chloride",
                "operator": ">=",
                "value": 0.5
            }
        })

    # Strategy is valid if:
    # 1. Final product has chloro-aryl
    # 2. At least one starting material has chloro-aryl
    # 3. At least 50% of molecule nodes have chloro-aryl
    result = (
        final_product_has_chloro_aryl
        and starting_material_has_chloro_aryl
        and preservation_ratio >= 0.5
    )

    return result, findings_json