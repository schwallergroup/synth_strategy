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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if a fluoro-aromatic group is preserved throughout the synthesis.

    It checks if a fluoro-aromatic group (a fluorine atom attached to any aromatic ring)
    is present in all non-starting materials in the main synthetic pathway. Starting
    materials (in_stock=True) are excluded from the check.
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

    # Track if all relevant molecules have fluorophenyl
    all_mols_have_fluorophenyl = True

    # Track if we found any non-starting materials
    found_non_starting_materials = False

    def has_fluorophenyl(mol_smiles):
        """Helper function to check if a molecule has a fluorophenyl group"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return False

        # Check for fluorine atoms attached to aromatic carbons
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "F":
                # Check if this fluorine is attached to an aromatic carbon
                neighbors = atom.GetNeighbors()
                if neighbors and neighbors[0].GetIsAromatic() and neighbors[0].GetSymbol() == "C":
                    return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal all_mols_have_fluorophenyl, found_non_starting_materials, findings_json

        if node["type"] == "mol":
            # Skip starting materials (in_stock=True)
            if not node.get("in_stock", False):
                mol_smiles = node.get("smiles", "")
                if mol_smiles:
                    found_non_materials_this_mol = False
                    if not found_non_starting_materials:
                        found_non_materials_this_mol = True
                    found_non_starting_materials = True

                    # Check if molecule has fluorophenyl group
                    if has_fluorophenyl(mol_smiles):
                        if "fluoro-aromatic group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("fluoro-aromatic group")
                    else:
                        print(f"Molecule without fluorophenyl group found: {mol_smiles}")
                        all_mols_have_fluorophenyl = False
                        # Add structural constraint for non-starting material without fluoro-aromatic group
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "non-starting_material_without_fluoro-aromatic_group",
                                "operator": "==",
                                "value": 0
                            }
                        })
                    if found_non_materials_this_mol:
                        # Add structural constraint for finding at least one non-starting material
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "non-starting_material",
                                "operator": ">=",
                                "value": 1
                            }
                        })

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    result = all_mols_have_fluorophenyl

    # If we didn't find any non-starting materials, return False
    if not found_non_starting_materials:
        print("No non-starting materials found in the synthesis route")
        result = False
        # If no non-starting materials, the constraint for >=1 non-starting material is not met
        # This is implicitly handled by the result being False, but we can explicitly add the constraint if needed.
        # For this specific case, the original function returns False, so we align with that.

    if result:
        print("Fluorophenyl group is preserved throughout the synthesis")

    return result, findings_json
