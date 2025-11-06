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
    This function detects a synthetic strategy involving late-stage nitro reduction to form an amine.
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

    has_nitro_reduction = False

    def dfs_traverse(node, current_depth=0):
        nonlocal has_nitro_reduction, findings_json

        # For reaction nodes, check if it's a nitro reduction
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Use depth from metadata if available, otherwise use calculated depth
            node_depth = node["metadata"].get("depth", current_depth)

            # Check if this is a late-stage reaction (within first 2 steps)
            if node_depth <= 2:
                print(f"Examining late-stage reaction at depth {node_depth}: {rsmi}")

                # Check if this is a nitro reduction reaction using the reaction checker
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Found nitro reduction reaction at depth {node_depth}")
                    has_nitro_reduction = True
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Reduction of nitro groups to amines",
                            "position": "late_stage",
                            "condition": "depth <= 2"
                        }
                    })

        # Recursively traverse children with new depth logic
        for child in node.get("children", []):
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, current_depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_nitro_reduction, findings_json
