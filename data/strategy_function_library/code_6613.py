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
    Detects a linear synthesis strategy that maintains a protected amine
    throughout the synthesis
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

    # Track linear paths with protected amines
    linear_paths_with_protected_amines = []

    def dfs_traverse(node, path=None, depth=0):
        nonlocal findings_json
        if path is None:
            path = []

        current_path = path.copy()
        current_path.append(node)

        # If this is a molecule node with no children, we've reached the end of a path
        if node["type"] == "mol" and not node.get("children"):
            # Check if this path has at least 2 reactions with protected amines
            protected_amine_reactions = 0
            protected_amine_present = False

            for i, node_in_path in enumerate(current_path):
                if node_in_path["type"] == "reaction":
                    rsmi = node_in_path["metadata"]["rsmi"]
                    product_str = rsmi.split(">")[-1]

                    # Check for protected amine in product using checker function
                    if checker.check_fg("Carbamic ester", product_str):
                        if "Carbamic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carbamic ester")
                        protected_amine_reactions += 1
                        protected_amine_present = True
                        print(f"Found protected amine in reaction at depth {depth-i}")
                    else:
                        # If we had a protected amine but lost it, reset the counter
                        if protected_amine_present:
                            protected_amine_present = False
                            protected_amine_reactions = 0

            # If we have at least 2 consecutive reactions with protected amines
            if protected_amine_reactions >= 2:
                linear_paths_with_protected_amines.append(current_path)
                if {"type": "count", "details": {"target": "reaction_product_has_Carbamic_ester", "operator": ">=", "value": 2, "scope": "consecutive"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reaction_product_has_Carbamic_ester", "operator": ">=", "value": 2, "scope": "consecutive"}})
                print(
                    f"Found linear path with {protected_amine_reactions} protected amine reactions"
                )

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not 'reaction' (e.g., 'mol'), increase depth
                new_depth = depth + 1
            # If current node is 'reaction', depth remains the same
            dfs_traverse(child, current_path, new_depth)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if at least one linear path with protected amines is found
    result = len(linear_paths_with_protected_amines) > 0
    return result, findings_json