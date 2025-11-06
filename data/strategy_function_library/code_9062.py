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
    Detects if a chlorophenyl group is maintained throughout the synthesis,
    indicating it's a key structural element that's preserved.
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

    # Store molecules with their depth in the synthesis route
    mol_data = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            # Only include non-starting materials
            if not node.get("in_stock", False):
                mol_data.append({"smiles": node["smiles"], "depth": depth})

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # If current node is not a reaction (e.g., chemical)
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Sort by depth (ascending) - final product first, early intermediates last
    mol_data.sort(key=lambda x: x["depth"])

    # Check if chlorophenyl is present in all significant molecules
    complex_mol_threshold = 0  # Include all molecules in the analysis

    complex_mols_with_chlorophenyl = 0
    total_complex_mols = 0

    for mol_info in mol_data:
        smiles = mol_info["smiles"]
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                total_complex_mols += 1
                # Check for chlorinated aromatic ring
                if checker.check_fg("Chlorobenzene", smiles):
                    complex_mols_with_chlorophenyl += 1
                    if "Chlorobenzene" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Chlorobenzene")
                    print(f"Found chlorophenyl in molecule: {smiles}")
        except Exception as e:
            print(f"Error processing molecule {smiles}: {e}")
            continue

    # Consider chlorophenyl maintained if at least 50% of molecules have it
    maintained = (
        total_complex_mols > 0 and complex_mols_with_chlorophenyl >= total_complex_mols * 0.5
    )

    if maintained:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "molecules_with_Chlorobenzene",
                "operator": ">=",
                "value": "0.5 * total_molecules"
            }
        })

    print(
        f"Chlorophenyl maintained: {maintained} ({complex_mols_with_chlorophenyl}/{total_complex_mols})"
    )
    return maintained, findings_json
