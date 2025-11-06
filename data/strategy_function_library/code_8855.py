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
    Detects if an aryl chloride group is present in a product and all of its non-commercially-available precursors.
    This strategy attempts to identify when an aryl chloride is maintained across a synthetic step.
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

    # Track if we've found an aryl chloride in the target molecule
    target_has_aryl_chloride = False
    # Track if any step loses the aryl chloride
    aryl_chloride_maintained = True

    def has_aromatic_chloride(mol_smiles):
        """Helper function to check if a molecule has an aromatic chloride"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if mol is None:
            return False

        # Verify that at least one of the aromatic halides is a chlorine
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "Cl" and atom.GetNeighbors()[0].GetIsAromatic():
                return True
        return False

    def dfs_traverse(node, depth=0, parent_has_aryl_cl=None):
        nonlocal target_has_aryl_chloride, aryl_chloride_maintained, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_aryl_chloride = has_aromatic_chloride(mol_smiles)

            # If this is the target molecule (depth 0)
            if depth == 0:
                target_has_aryl_chloride = has_aryl_chloride
                print(f"Target molecule has aryl chloride: {target_has_aryl_chloride}")
                if target_has_aryl_chloride:
                    findings_json["atomic_checks"]["functional_groups"].append("aryl chloride")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "aryl chloride",
                            "position": "last_stage"
                        }
                    })

                # If target doesn't have aryl chloride, no need to check maintenance
                if not target_has_aryl_chloride:
                    return

            # For non-target molecules, check if aryl chloride is maintained
            elif parent_has_aryl_cl is not None:
                # Only check non-starting materials
                if parent_has_aryl_cl and not has_aryl_chloride and not node.get("in_stock", False):
                    print(f"Aryl chloride lost at depth {depth} in molecule: {mol_smiles}")
                    aryl_chloride_maintained = False
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "destruction of aryl chloride"
                        }
                    })

            # Pass down whether this molecule has aryl chloride
            for child in node.get("children", []):
                # New logic: depth increases only from chemical to reaction
                new_depth = depth + 1 if node["type"] != "reaction" else depth
                dfs_traverse(child, new_depth, has_aryl_chloride)

        elif node["type"] == "reaction":
            # For reaction nodes, we just pass through the parent's aryl chloride status
            for child in node.get("children", []):
                # New logic: depth increases only from chemical to reaction
                new_depth = depth + 1 if node["type"] != "reaction" else depth
                dfs_traverse(child, new_depth, parent_has_aryl_cl)

    # Start traversal
    dfs_traverse(route)

    result = False
    # If target doesn't have aryl chloride, return False (no maintenance to check)
    if not target_has_aryl_chloride:
        print("Target molecule doesn't have aryl chloride, so maintenance check is not applicable")
        result = False
    else:
        print(f"Aryl chloride maintained throughout synthesis: {aryl_chloride_maintained}")
        result = aryl_chloride_maintained

    return result, findings_json