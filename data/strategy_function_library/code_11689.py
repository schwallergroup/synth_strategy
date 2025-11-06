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
    This function detects if the synthetic route maintains a polycyclic scaffold throughout
    a linear synthesis. It checks for the presence of multiple rings in intermediates.
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

    # Track molecules and their properties at each depth
    molecules_by_depth = {}

    # Track linearity and polycyclic scaffold presence
    step_count = 0
    polycyclic_steps = 0

    # Function to check if a molecule has a polycyclic scaffold
    def has_polycyclic_scaffold(mol_smiles):
        try:
            mol = Chem.MolFromSmiles(mol_smiles)
            if not mol:
                return False

            # Get ring information
            ring_info = mol.GetRingInfo()

            # Get all rings
            rings = ring_info.AtomRings()
            if len(rings) < 2:  # Changed from 3 to 2 to include bicyclic systems
                return False

            # Check if rings are connected (share atoms)
            # Build a graph of ring connections
            ring_graph = {}
            for i in range(len(rings)):
                ring_graph[i] = set()
                for j in range(len(rings)):
                    if i != j:
                        # Check if rings i and j share any atoms
                        if set(rings[i]).intersection(set(rings[j])):
                            ring_graph[i].add(j)

            # Check if we have at least 2 connected rings (changed from 3)
            visited = set()

            def dfs_rings(ring_idx):
                visited.add(ring_idx)
                for neighbor in ring_graph[ring_idx]:
                    if neighbor not in visited:
                        dfs_rings(neighbor)

            # Start DFS from first ring
            if rings:
                dfs_rings(0)

            # If we can visit at least 2 rings, we have a polycyclic system (changed from 3)
            if len(visited) >= 2:
                # Add a generic polycyclic ring system finding if detected
                if "Polycyclic scaffold" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("Polycyclic scaffold")
                return True
            return False

        except Exception as e:
            print(f"Error checking polycyclic scaffold: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal step_count, polycyclic_steps, findings_json

        if node["type"] == "mol":
            # Store molecule at this depth
            if depth not in molecules_by_depth:
                molecules_by_depth[depth] = []
            molecules_by_depth[depth].append(node["smiles"])

            # Check if this molecule has a polycyclic scaffold
            if has_polycyclic_scaffold(node["smiles"]):
                print(f"Found polycyclic scaffold in molecule at depth {depth}: {node['smiles']}")

        elif node["type"] == "reaction":
            step_count += 1

            # Check if product has polycyclic scaffold
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                product = rsmi.split(">")[-1]

                if has_polycyclic_scaffold(product):
                    polycyclic_steps += 1
                    print(f"Found polycyclic scaffold in reaction product at depth {depth}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from root
    dfs_traverse(route)

    # Check linearity by examining molecule counts at each depth
    is_linear = True
    for depth, molecules in molecules_by_depth.items():
        # In a linear synthesis, we should have at most 3 molecules at any depth
        # (allowing for more branching)
        if len(molecules) > 3:  # Changed from 2 to 3
            is_linear = False
            print(f"Non-linear synthesis detected at depth {depth} with {len(molecules)} molecules")
            # Add structural constraint finding
            if {"type": "count", "details": {"target": "molecules_per_step", "operator": "<=", "value": 3}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "count", "details": {"target": "molecules_per_step", "operator": "<=", "value": 3}})
            break

    # Count total number of leaf nodes (starting materials)
    leaf_count = 0

    def count_leaves(node):
        nonlocal leaf_count
        if node["type"] == "mol" and not node.get("children", []):
            leaf_count += 1
        for child in node.get("children", []):
            count_leaves(child)

    count_leaves(route)

    # A linear synthesis can have more starting materials in practice
    if leaf_count > 6: # Changed from 3 to 6
        is_linear = False
        # Add structural constraint finding
        if {"type": "count", "details": {"target": "leaf_nodes", "operator": "<=", "value": 6}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "leaf_nodes", "operator": "<=", "value": 6}})

    print(f"Total steps: {step_count}, Polycyclic steps: {polycyclic_steps}")
    print(f"Is linear: {is_linear}, Leaf count: {leaf_count}")

    # Check if most steps maintain a polycyclic scaffold (at least 70% of steps)
    has_polycyclic_scaffold_throughout = step_count > 0 and polycyclic_steps / step_count >= 0.7
    if has_polycyclic_scaffold_throughout:
        if {"type": "count", "details": {"target": "ratio_of_steps_with_polycyclic_product", "operator": ">=", "value": 0.7}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "ratio_of_steps_with_polycyclic_product", "operator": ">=", "value": 0.7}})

    # Check if the target molecule (depth 0) has a polycyclic scaffold
    target_has_polycyclic = False
    if 0 in molecules_by_depth and molecules_by_depth[0]:
        target_has_polycyclic = has_polycyclic_scaffold(molecules_by_depth[0][0])
        if target_has_polycyclic:
            if {"type": "positional", "details": {"target": "polycyclic_scaffold", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "polycyclic_scaffold", "position": "last_stage"}})

    # For a true linear synthesis with polycyclic scaffold:
    # 1. The synthesis must be linear
    # 2. The target molecule must have a polycyclic scaffold
    # 3. At least 70% of steps must maintain the polycyclic scaffold
    result = is_linear and target_has_polycyclic and has_polycyclic_scaffold_throughout
    print(f"Final result: {result}")
    print(f"- Linear synthesis: {is_linear}")
    print(f"- Target has polycyclic scaffold: {target_has_polycyclic}")
    print(f"- Maintains polycyclic scaffold throughout: {has_polycyclic_scaffold_throughout}")

    return result, findings_json
