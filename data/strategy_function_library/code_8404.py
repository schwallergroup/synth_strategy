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
    Detects if fluorophenyl groups are maintained throughout the synthesis.
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

    depths_with_fluorophenyl = set()
    mol_depths = set()
    max_depth = 0
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal depths_with_fluorophenyl, mol_depths, max_depth, findings_json

        if node["type"] == "mol":
            mol_depths.add(depth)
            mol_smiles = node["smiles"]

            mol = Chem.MolFromSmiles(mol_smiles)
            if mol:
                # Look for aromatic carbon connected to fluorine
                fluorophenyl_pattern = Chem.MolFromSmarts("c-F")
                if mol.HasSubstructMatch(fluorophenyl_pattern):
                    depths_with_fluorophenyl.add(depth)
                    if "fluorophenyl" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("fluorophenyl")
                    print(
                        f"Fluorophenyl group detected at depth {depth} in molecule: {mol_smiles}"
                    )

        max_depth = max(max_depth, depth)

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not 'reaction' (e.g., 'mol'), increase depth
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)

    # Check if fluorophenyl groups are present at all molecule depths
    if not mol_depths:
        print("No molecule nodes found in the route")
        return False, findings_json

    all_mol_depths_with_fluorophenyl = all(
        depth in depths_with_fluorophenyl for depth in mol_depths
    )

    if all_mol_depths_with_fluorophenyl:
        print("Fluorophenyl groups maintained throughout the synthesis")
        result = True
    else:
        missing_depths = [d for d in mol_depths if d not in depths_with_fluorophenyl]
        print(f"Fluorophenyl groups missing at depths: {missing_depths}")
        # If fluorophenyl groups are NOT maintained, it means the negation constraint is met.
        # The original strategy JSON defines a 'negation' constraint for 'absence_of_fluorophenyl_in_any_molecule'.
        # If 'all_mol_depths_with_fluorophenyl' is False, it implies absence at some depth.
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "absence_of_fluorophenyl_in_any_molecule"
            }
        })

    return result, findings_json