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


CC_COUPLING_REACTIONS = [
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl OTf",
    "Sonogashira alkyne_aryl OTf",
    "Sonogashira acetylene_alkenyl halide",
    "Sonogashira alkyne_alkenyl halide",
    "Sonogashira acetylene_alkenyl OTf",
    "Sonogashira alkyne_alkenyl OTf",
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Negishi coupling",
    "Stille reaction_aryl",
    "Stille reaction_vinyl",
    "Heck terminal vinyl",
    "Kumada cross-coupling",
    "Hiyama-Denmark Coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects an early-stage C-C bond formation strategy. The function checks if any of the first three steps of a synthesis is a cross-coupling reaction defined in the CC_COUPLING_REACTIONS list.
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

    found_cc_coupling = False

    # Pre-traversal to find the maximum depth of the synthesis route
    max_depth = 0
    if route:
        q = [(route, 0)]
        visited_nodes = {id(route)}
        while q:
            node, depth = q.pop(0)
            if node.get("type") == "reaction":
                max_depth = max(max_depth, depth)
            for child in node.get("children", []):
                if id(child) not in visited_nodes:
                    q.append((child, depth + 1))
                    visited_nodes.add(id(child))

    def dfs_traverse(node, depth, max_depth):
        nonlocal found_cc_coupling, findings_json

        # "Early-stage" is defined as the first 3 steps (steps at depth >= max_depth - 2).
        # This check is robust for syntheses of any length.
        if node["type"] == "reaction" and depth >= max_depth - 2:
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is a C-C coupling reaction from the predefined list
                for reaction_type in CC_COUPLING_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        found_cc_coupling = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        # Add the structural constraint if a CC coupling is found in the early stages
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "any_cc_coupling_reaction",
                                "position": "first_three_stages"
                            }
                        })
                        break

        # Traverse children
        for child in node.get("children", []):
            if not found_cc_coupling:
                # New logic for depth calculation
                new_depth = depth
                if node["type"] != "reaction": # If current node is chemical, depth increases
                    new_depth = depth + 1
                # If current node is reaction, depth remains the same
                dfs_traverse(child, new_depth, max_depth)

    # Start traversal from the root
    if route:
        dfs_traverse(route, 0, max_depth)

    return found_cc_coupling, findings_json
