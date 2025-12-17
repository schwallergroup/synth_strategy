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
    This function detects a synthetic strategy that maintains a furan ring throughout the synthesis.
    It checks if the target molecule and all key intermediates (non-starting materials) contain a furan ring.
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

    # Track if the main synthetic pathway preserves the furan ring
    main_pathway_preserves_furan = True
    # Count molecules in the main synthetic pathway (excluding starting materials)
    main_pathway_mol_count = 0

    def dfs_traverse(node, is_main_pathway=True, depth=0):
        nonlocal main_pathway_preserves_furan, main_pathway_mol_count, findings_json

        if node["type"] == "mol" and "smiles" in node:
            # Skip checking starting materials
            is_starting_material = node.get("in_stock", False)

            # Only check molecules in the main pathway that aren't starting materials
            if is_main_pathway and not is_starting_material:
                main_pathway_mol_count += 1
                has_furan = checker.check_ring("furan", node["smiles"])
                if has_furan:
                    # If furan is found, record it. The strategy is about *maintaining* it.
                    # We record it here to show that the check for furan was performed and passed for this molecule.
                    if "furan" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("furan")
                else:
                    main_pathway_preserves_furan = False

        # Determine the new depth for the recursive call
        new_depth = depth
        if node["type"] != "reaction":  # If current node is chemical, depth increases
            new_depth = depth + 1

        # Traverse children
        for i, child in enumerate(node.get("children", [])):
            # Only the first child in a reaction node is considered part of the main pathway
            # This assumes the first child is the primary reactant in retrosynthetic direction
            if node["type"] == "reaction":
                is_child_main_pathway = (i == 0) and is_main_pathway
            else:
                is_child_main_pathway = is_main_pathway

            dfs_traverse(child, is_child_main_pathway, new_depth)

    # Start traversal
    dfs_traverse(route)

    result = main_pathway_preserves_furan and main_pathway_mol_count > 0

    # Record structural constraints based on the final result
    if result:
        # If the main pathway preserves furan, it means the negation constraint is met (no intermediate lacked furan)
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "main_pathway_intermediate_lacks_furan"
            }
        })
        # If main_pathway_mol_count > 0, it means the count constraint is met
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "main_pathway_intermediate",
                "operator": ">",
                "value": 0
            }
        })

    return result, findings_json
