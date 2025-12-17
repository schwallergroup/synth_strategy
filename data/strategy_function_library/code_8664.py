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
    Detects if a fluorinated aromatic system is maintained throughout the synthesis.
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

    # Track if we've found a fluorinated aromatic system at each depth
    fluoro_aromatic_at_depth = {}
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            print(f"Processing molecule at depth {depth}: {mol_smiles}")

            # Check for fluorinated aromatic
            has_fluoro_aromatic = checker.check_fg("Aromatic halide", mol_smiles)

            # Verify it's specifically a fluorinated aromatic
            if has_fluoro_aromatic:
                mol = Chem.MolFromSmiles(mol_smiles)
                is_actual_fluoro_aromatic = False
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() == "F" and any(
                        neigh.GetIsAromatic() for neigh in atom.GetNeighbors()
                    ):
                        is_actual_fluoro_aromatic = True
                        break
                if is_actual_fluoro_aromatic:
                    if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                has_fluoro_aromatic = is_actual_fluoro_aromatic

            # Record if this depth has a fluorinated aromatic
            if depth not in fluoro_aromatic_at_depth:
                fluoro_aromatic_at_depth[depth] = has_fluoro_aromatic
            else:
                fluoro_aromatic_at_depth[depth] = (
                    fluoro_aromatic_at_depth[depth] or has_fluoro_aromatic
                )

            print(f"Depth {depth}: {'Has' if has_fluoro_aromatic else 'No'} fluorinated aromatic")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a 'mol' node
            next_depth = depth + 1

        # Traverse children (reactants)
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if fluorinated aromatic is present at all depths
    if not fluoro_aromatic_at_depth:
        print("No molecules found in route")
        result = False
        return result, findings_json

    # Get the maximum depth
    max_depth = max(fluoro_aromatic_at_depth.keys())

    # Check if fluorinated aromatic exists at each depth from 0 to max_depth
    all_depths_have_fluoro_aromatic = True
    for d in range(max_depth + 1):
        if d not in fluoro_aromatic_at_depth or not fluoro_aromatic_at_depth[d]:
            print(f"No fluorinated aromatic at depth {d}")
            all_depths_have_fluoro_aromatic = False
            break

    if all_depths_have_fluoro_aromatic:
        print(f"Fluorinated aromatic maintained through all {max_depth+1} depths")
        result = True
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Aromatic halide",
                "position": "all_stages"
            }
        })
    else:
        result = False

    return result, findings_json