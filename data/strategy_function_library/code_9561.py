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
    Detects if the synthesis involves a trifluoromethyl-containing aromatic fragment
    that persists throughout the synthesis.
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

    result = False

    # Track depths where we find trifluoromethyl groups and their atom mappings
    depths_with_cf3 = {}
    cf3_atom_maps = {}

    def has_aromatic_cf3(mol_smiles):
        """Check if the molecule has a CF3 group attached to an aromatic ring"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return False

        # Find all CF3 groups
        cf3_atoms = mol.GetSubstructMatches(Chem.MolFromSmarts("C(F)(F)F"))

        # Check if any CF3 is attached to an aromatic carbon
        for cf3_match in cf3_atoms:
            cf3_carbon = cf3_match[0]
            for neighbor in mol.GetAtomWithIdx(cf3_carbon).GetNeighbors():
                if neighbor.GetIsAromatic():
                    return True
        return False

    def get_cf3_atom_maps(mol_smiles):
        """Extract atom mapping numbers for CF3 groups attached to aromatic rings"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return []

        result = []
        # Find all CF3 groups
        cf3_atoms = mol.GetSubstructMatches(Chem.MolFromSmarts("C(F)(F)F"))

        for cf3_match in cf3_atoms:
            cf3_carbon = cf3_match[0]
            for neighbor in mol.GetAtomWithIdx(cf3_carbon).GetNeighbors():
                if neighbor.GetIsAromatic():
                    # Get atom mapping of the CF3 carbon
                    atom_map = (
                        str(mol.GetAtomWithIdx(cf3_carbon).GetAtomMapNum())
                        if mol.GetAtomWithIdx(cf3_carbon).GetAtomMapNum() > 0
                        else None
                    )
                    if atom_map:
                        result.append(atom_map)
        return result

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for trifluoromethyl group attached to aromatic ring
            if has_aromatic_cf3(mol_smiles):
                print(f"Found aromatic trifluoromethyl at depth {depth} in molecule: {mol_smiles}")
                if "trifluoromethyl group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("trifluoromethyl group")

                # Get atom mappings for CF3 groups
                atom_maps = get_cf3_atom_maps(mol_smiles)
                print(f"Atom maps for CF3 at depth {depth}: {atom_maps}")

                for atom_map in atom_maps:
                    if atom_map not in cf3_atom_maps:
                        cf3_atom_maps[atom_map] = set()
                    cf3_atom_maps[atom_map].add(depth)

                # Also track depths
                if depth not in depths_with_cf3:
                    depths_with_cf3[depth] = 0
                depths_with_cf3[depth] += 1

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    print(f"CF3 atom maps and their depths: {cf3_atom_maps}")
    print(f"Depths with CF3: {depths_with_cf3}")

    # Check if any CF3 group persists through at least 3 steps
    for atom_map, depths in cf3_atom_maps.items():
        if len(depths) >= 3:
            print(f"CF3 with atom map {atom_map} persists through {len(depths)} steps")
            result = True
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "persistence_of_single_aromatic_trifluoromethyl_group", "operator": ">=", "value": 3}})
            # No need to break here, as we want to capture all relevant findings

    # If we have CF3 groups at 3 or more depths and they span a significant portion of the synthesis
    if len(depths_with_cf3) >= 3:
        depths_list = sorted(depths_with_cf3.keys())

        # Check if CF3 groups span a significant portion of the synthesis
        synthesis_span = max(depths_list) - min(depths_list)
        if synthesis_span >= 2:
            print(
                f"Found CF3 groups spanning {synthesis_span+1} depths from {min(depths_list)} to {max(depths_list)}"
            )

    # Even without atom mapping, if we have CF3 at 4 or more different depths, it's likely persistent
    if len(depths_with_cf3) >= 4:
        print(
            f"Found CF3 groups at {len(depths_with_cf3)} different depths, indicating persistence"
        )
        result = True
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "synthesis_steps_with_aromatic_trifluoromethyl_group", "operator": ">=", "value": 4}})

    return result, findings_json