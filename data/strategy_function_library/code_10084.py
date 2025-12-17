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
    Detects if a fluorinated aromatic moiety is maintained throughout the synthesis.
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

    def has_fluorinated_aromatic(mol_smiles):
        """Helper function to check for fluorinated aromatic"""
        if not mol_smiles:
            return False

        # Check for aromatic halide
        if not checker.check_fg("Aromatic halide", mol_smiles):
            return False
        else:
            if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

        # Verify it contains fluorine
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return False

        # Check specifically for aromatic carbon with fluorine attached
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic() and atom.GetSymbol() == "C":
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == "F":
                        return True
        return False

    # First check if the target molecule has a fluorinated aromatic
    if route["type"] == "mol":
        target_smiles = route.get("smiles", "")
        if not target_smiles:
            print("Target molecule has no SMILES")
            return False, findings_json

        if not has_fluorinated_aromatic(target_smiles):
            print("Target molecule does not have a fluorinated aromatic")
            return False, findings_json
        else:
            # If target has it, it means 'Aromatic halide' is present at the last stage
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "Aromatic halide",
                    "position": "last_stage"
                }
            })

    # Track paths where fluorinated aromatic is maintained
    valid_paths = []

    def dfs_traverse(node, path=None, has_f_aromatic=False):
        """
        Traverse the synthesis route and track if fluorinated aromatic is maintained.
        In retrosynthetic direction: from target to starting materials.
        """
        nonlocal findings_json # Declare findings_json as nonlocal to modify it

        if path is None:
            path = []

        current_path = path + [node]

        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if not smiles:
                return False

            current_has_f_aromatic = has_fluorinated_aromatic(smiles)

            # If we're tracking a fluorinated aromatic and it disappears, this path is invalid
            if has_f_aromatic and not current_has_f_aromatic:
                # This indicates a loss, so we don't add the negation constraint here
                # The overall result will be False, which implicitly means the negation condition is met.
                return False

            # If this is a starting material (leaf node)
            if node.get("in_stock", False) or not node.get("children", []):
                # If we've been tracking a fluorinated aromatic and it's still present, this is a valid path
                if has_f_aromatic:
                    valid_paths.append(current_path)
                    return True
                return False

        elif node["type"] == "reaction":
            # The logic in this block was redundant. The traversal to children handles
            # the check correctly when it reaches the next molecule node. The reaction
            # node simply acts as a connector in the graph.
            pass

        # Continue traversal to children
        if not node.get("children", []):
            return False

        # If any child branch maintains the fluorinated aromatic, this branch is valid
        return any(
            dfs_traverse(child, current_path, has_f_aromatic) for child in node.get("children", [])
        )

    # Start traversal from the root with fluorinated aromatic tracking enabled
    # (since we verified the target has it)
    result = dfs_traverse(route, None, True)

    if result:
        print(
            f"Detected maintenance of fluorinated aromatic throughout synthesis ({len(valid_paths)} valid paths)"
        )
        # If result is True, it means there was NO loss of Aromatic halide during synthesis
        # So, the negation constraint is NOT met (which is good for the strategy)
        # We only add structural constraints if they are DETECTED.
        # The strategy implies we want to detect if the negation condition (loss) is NOT met.
        # So, if result is True, it means 'loss_of_Aromatic_halide_during_synthesis' did NOT happen.
        # The current JSON structure for structural_constraints is for *detected* constraints.
        # If the strategy is 'negation' and the condition is met (i.e., the loss did NOT occur),
        # then we should record that the negation constraint was satisfied.
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "loss_of_Aromatic_halide_during_synthesis"
            }
        })

    else:
        print("No path maintains fluorinated aromatic throughout synthesis")
        # If result is False, it means there WAS a loss of Aromatic halide during synthesis.
        # In this case, the negation constraint (that there should be no loss) is NOT satisfied.
        # We don't add it to findings_json because it wasn't 'found' in the positive sense.

    return result, findings_json
