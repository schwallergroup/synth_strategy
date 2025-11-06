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
    This function detects if key functional groups (aniline and nitrile) are preserved
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

    # Track presence of functional groups at each step
    steps_with_aniline = set()
    steps_with_nitrile = set()
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, findings_json

        if node["type"] == "mol":
            # Update max depth for molecules
            max_depth = max(max_depth, depth)

            # Check for functional groups in molecule
            mol_smiles = node["smiles"]

            # Check for aniline
            if checker.check_fg("Aniline", mol_smiles):
                steps_with_aniline.add(depth)
                if "Aniline" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Aniline")
                print(f"Detected aniline at depth {depth} in molecule: {mol_smiles}")

            # Check for nitrile
            if checker.check_fg("Nitrile", mol_smiles):
                steps_with_nitrile.add(depth)
                if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                print(f"Detected nitrile at depth {depth} in molecule: {mol_smiles}")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only from chemical to reaction node
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node['type'] is 'mol' or 'chemical'
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    print(f"Max depth: {max_depth}")
    print(f"Steps with aniline: {sorted(steps_with_aniline)}")
    print(f"Steps with nitrile: {sorted(steps_with_nitrile)}")

    # Check if functional groups are preserved from start to end
    aniline_preserved = 0 in steps_with_aniline and max_depth in steps_with_aniline
    nitrile_preserved = 0 in steps_with_nitrile and max_depth in steps_with_nitrile

    strategy_present = (
        aniline_preserved or nitrile_preserved
    )

    if strategy_present:
        if aniline_preserved:
            print("Detected aniline preservation strategy")
            # Add structural constraints for Aniline preservation
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Aniline", "position": "first_stage"}})
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Aniline", "position": "last_stage"}})
        if nitrile_preserved:
            print("Detected nitrile preservation strategy")
            # Add structural constraints for Nitrile preservation
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Nitrile", "position": "first_stage"}})
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Nitrile", "position": "last_stage"}})

    return strategy_present, findings_json
