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
    Detects a synthesis pathway of at least three steps where a nitro group is present in every molecule along the main synthetic path.
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

    # Track the main synthetic pathway
    main_pathway = []
    result = False

    def dfs_traverse(node, depth=0, path=None):
        nonlocal main_pathway, result, findings_json
        """Traverse the synthesis tree and identify molecules with nitro groups"""
        if path is None:
            path = []

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if it has a nitro group
            has_nitro = checker.check_fg("Nitro group", mol_smiles)

            # Add to path if it has a nitro group and is not in stock
            if has_nitro:
                if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

                # Create a copy of the current path and add this molecule
                current_path = path + [
                    {"smiles": mol_smiles, "depth": depth, "in_stock": node.get("in_stock", False)}
                ]

                # If this is a leaf node (no children) or in_stock, consider this a complete path
                if not node.get("children", []) or node.get("in_stock", False):
                    # Only consider paths with non-stock molecules (excluding the starting material)
                    non_stock_mols = [m for m in current_path if not m["in_stock"]]
                    if (
                        len(non_stock_mols) >= 3
                    ):  # Need at least 3 molecules for a meaningful synthesis
                        # If this path is longer than our current best, update it
                        if len(non_stock_mols) > len(
                            [m for m in main_pathway if not m["in_stock"]]
                        ):
                            main_pathway.clear()
                            main_pathway.extend(current_path)

                # Continue traversing children
                for child in node.get("children", []):
                    # Depth increases when going from chemical to reaction
                    dfs_traverse(child, depth + 1, current_path)
            else:
                # If this molecule doesn't have a nitro group, don't include it in the path
                # but still traverse its children (might be a branch that doesn't preserve nitro)
                for child in node.get("children", []):
                    # Depth increases when going from chemical to reaction
                    dfs_traverse(child, depth + 1, path)

        elif node["type"] == "reaction":
            # For reaction nodes, just continue traversing
            for child in node.get("children", []):
                # Depth remains the same when going from reaction to chemical
                dfs_traverse(child, depth, path)

    # Start traversal from the root node
    dfs_traverse(route)

    # If we didn't find any pathway with at least 3 molecules, return False
    if len(main_pathway) < 3:
        result = False
    else:
        # Sort molecules by depth (from target to starting materials)
        main_pathway.sort(key=lambda x: x["depth"])

        # Filter out in_stock molecules except the last one (starting material)
        filtered_pathway = []
        for i, mol_info in enumerate(main_pathway):
            if not mol_info["in_stock"] or i == len(main_pathway) - 1:
                filtered_pathway.append(mol_info)

        # Need at least 3 molecules for a meaningful synthesis
        if len(filtered_pathway) >= 3:
            result = True
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "synthesis_steps",
                    "operator": ">=",
                    "value": 3
                }
            })
            # All molecules on the main path must have a nitro group, which is implicitly handled by how main_pathway is built.
            # If result is True, it means all molecules in the filtered_pathway had a nitro group.
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "molecules_without_nitro_group_on_main_path",
                    "operator": "==",
                    "value": 0
                }
            })
        else:
            result = False

    return result, findings_json
