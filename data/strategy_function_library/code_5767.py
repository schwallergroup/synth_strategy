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
    Detects if a trifluoromethyl ether group is present in a high proportion (>= 50%) of molecules throughout a synthesis, as a heuristic for a potential preservation strategy.
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

    steps_with_trifluoromethyl_ether = 0
    total_mol_steps = 0

    # Track molecules with the functional group by their depth
    depths_with_trifluoromethyl = []

    def dfs_traverse(node, depth=0):
        nonlocal steps_with_trifluoromethyl_ether, total_mol_steps, findings_json

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            total_mol_steps += 1

            # Use the checker function to detect trifluoromethyl ether group
            if checker.check_fg("Trifluoromethyl ether", mol_smiles):
                steps_with_trifluoromethyl_ether += 1
                depths_with_trifluoromethyl.append(depth)
                findings_json["atomic_checks"]["functional_groups"].append("Trifluoromethyl ether")
                print(f"Trifluoromethyl ether found at depth {depth}, SMILES: {mol_smiles}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a 'chemical' node or 'mol'
            next_depth = depth + 1

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Calculate percentage of molecules with trifluoromethyl ether
    preservation_percentage = 0
    if total_mol_steps > 0:
        preservation_percentage = (steps_with_trifluoromethyl_ether / total_mol_steps) * 100

    print(f"Total molecules: {total_mol_steps}")
    print(f"Molecules with trifluoromethyl ether: {steps_with_trifluoromethyl_ether}")
    print(f"Preservation percentage: {preservation_percentage:.2f}%")
    print(f"Depths with trifluoromethyl ether: {depths_with_trifluoromethyl}")

    result = False
    # Check if trifluoromethyl ether is present in at least 50% of molecules
    # and appears at multiple depths (indicating preservation)
    if (
        total_mol_steps > 0
        and steps_with_trifluoromethyl_ether / total_mol_steps >= 0.5
        and len(depths_with_trifluoromethyl) >= 2
    ):
        print("Trifluoromethyl ether preservation strategy detected")
        result = True
        # Add structural constraints to findings_json
        if total_mol_steps > 0 and steps_with_trifluoromethyl_ether / total_mol_steps >= 0.5:
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "proportion_of_molecules_with_Trifluoromethyl_ether",
                    "operator": ">=",
                    "value": 0.5
                }
            })
        if len(depths_with_trifluoromethyl) >= 2:
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "Trifluoromethyl ether",
                    "operator": ">=",
                    "value": 2
                }
            })

    if not result:
        print("No trifluoromethyl ether preservation strategy detected")
    
    return result, findings_json
