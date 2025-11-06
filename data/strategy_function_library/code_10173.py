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


HETEROCYCLES_OF_INTEREST = ["thiophene", "furan", "pyrrole", "pyridine", "imidazole", "oxazole", "thiazole"]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks if a specific heterocycle, identified in the final product, is also present in the starting material and at least 70% of all molecular intermediates. The heterocycle must be one of the types defined in the HETEROCYCLES_OF_INTEREST list.
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

    # Track depths and types of nodes
    depths_with_heterocycle = set()
    mol_depths = set()
    max_depth = 0

    # Track which heterocycle is being maintained
    maintained_heterocycle = None

    def dfs_traverse(node, depth=0):
        nonlocal depths_with_heterocycle, mol_depths, max_depth, maintained_heterocycle, findings_json

        if node["type"] == "mol":
            mol_depths.add(depth)
            max_depth = max(max_depth, depth)

            # Check for heterocycles in the molecule
            for heterocycle in HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(heterocycle, node["smiles"]):
                    findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                    # If this is the first heterocycle we find in the final product, set it as the one to maintain
                    if maintained_heterocycle is None and depth == 0:
                        maintained_heterocycle = heterocycle
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "any_heterocycle_of_interest", "position": "last_stage"}})

                    # If we've already found a heterocycle to maintain, only count this one if it matches
                    if maintained_heterocycle == heterocycle:
                        depths_with_heterocycle.add(depth)

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not 'reaction' (e.g., 'chemical' or 'mol')
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # No heterocycle found in the final product to track
    if maintained_heterocycle is None:
        return False, findings_json

    # Check if the heterocycle is present at the end (depth=0) and start (depth=max_depth) of synthesis
    end_has_heterocycle = 0 in depths_with_heterocycle
    start_has_heterocycle = max_depth in depths_with_heterocycle and max_depth in mol_depths

    if start_has_heterocycle:
        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "maintained_heterocycle", "position": "first_stage"}})

    # Calculate the percentage of molecule nodes that contain the heterocycle
    mol_depths_with_heterocycle = depths_with_heterocycle.intersection(mol_depths)
    continuity_percentage = len(mol_depths_with_heterocycle) / len(mol_depths) if mol_depths else 0

    if continuity_percentage >= 0.7:
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "molecules_with_maintained_heterocycle", "operator": ">=", "value": 0.7, "is_ratio": True}})

    if max_depth > 0:
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "synthesis_steps", "operator": ">", "value": 0}})

    # We consider the heterocycle maintained if:
    # 1. It's present at the start and end of synthesis
    # 2. It's present in at least 70% of molecule nodes
    # 3. The synthesis has at least one step (max_depth > 0)
    is_maintained = (
        start_has_heterocycle
        and end_has_heterocycle
        and continuity_percentage >= 0.7
        and max_depth > 0
    )

    return is_maintained, findings_json
