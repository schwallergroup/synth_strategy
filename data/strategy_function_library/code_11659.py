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
    Detects if a diaryl ether motif is maintained throughout the synthesis.
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

    # Track if we found the pattern at each step
    steps_with_diaryl_ether = 0
    total_steps = 0

    def has_diaryl_ether(mol_smiles):
        """Helper function to check if a molecule contains a diaryl ether motif"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return False

        # A specific SMARTS pattern is the most robust way to find a diaryl ether.
        # 'c-O-c' finds an oxygen connected to two aromatic carbons.
        diaryl_ether_pattern = Chem.MolFromSmarts("c-O-c")
        if not diaryl_ether_pattern:
            return False

        matches = mol.GetSubstructMatches(diaryl_ether_pattern)

        if matches:
            # Verify that the ether links two distinct aromatic rings,
            # which excludes fused systems like dibenzofuran.
            ring_info = mol.GetRingInfo()
            for match in matches:
                c1_idx, o_idx, c2_idx = match
                # Check if the two aromatic carbons are not in the same ring system.
                if not ring_info.AreAtomsInSameRing(c1_idx, c2_idx):
                    return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal steps_with_diaryl_ether, total_steps, findings_json

        if node["type"] == "mol" and not node.get("in_stock", False):
            total_steps += 1
            mol_smiles = node["smiles"]

            if has_diaryl_ether(mol_smiles):
                steps_with_diaryl_ether += 1
                # Record the finding of diaryl ether
                if "diaryl ether" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("diaryl ether")

        # Determine the new depth for children based on the current node's type
        new_depth = depth
        if node["type"] != "reaction":  # If current node is chemical, depth increases
            new_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Return True if diaryl ether is present in all non-starting material molecules except at most one
    result = total_steps > 0 and steps_with_diaryl_ether >= total_steps - 1

    # Record structural constraint if the condition is met
    if total_steps > 0 and (total_steps - steps_with_diaryl_ether) <= 1:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "absence of diaryl ether in intermediate",
                "operator": "<=",
                "value": 1
            }
        })

    return result, findings_json
