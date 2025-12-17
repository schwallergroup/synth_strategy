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
    Detects if the synthesis uses fluorinated aromatic building blocks.
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

    has_fluorinated_aromatics = False

    def dfs_traverse(node, depth=0):
        nonlocal has_fluorinated_aromatics, findings_json

        if node["type"] == "mol":
            # Check molecules for fluorinated aromatics
            mol_smiles = node["smiles"]
            mol = Chem.MolFromSmiles(mol_smiles)

            if mol is None:
                print(f"Warning: Could not parse molecule SMILES: {mol_smiles}")
                return

            # Check if the molecule has aromatic atoms
            has_aromatic = False
            for atom in mol.GetAtoms():
                if atom.GetIsAromatic():
                    has_aromatic = True
                    break

            if has_aromatic:
                # Check for aromatic halides (specifically fluorine)
                if checker.check_fg("Aromatic halide", mol_smiles) and mol.HasSubstructMatch(
                    Chem.MolFromSmarts("c-F")
                ):
                    has_fluorinated_aromatics = True
                    findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                    # Note: "fluorinated aromatic halide building block" is a conceptual target for structural constraint, not an atomic check name.
                    print(
                        f"Detected fluorinated aromatic building block (aromatic halide): {mol_smiles}"
                    )

                # Check for trifluoro groups on aromatics
                trifluoro_aromatic = Chem.MolFromSmarts("c-[#6](-F)(-F)-F")
                if mol.HasSubstructMatch(trifluoro_aromatic):
                    has_fluorinated_aromatics = True
                    findings_json["atomic_checks"]["functional_groups"].append("aromatic trifluoromethyl") # Assuming this maps to 'aromatic trifluoromethyl' in the strategy
                    # Note: "aromatic trifluoromethyl building block" is a conceptual target for structural constraint, not an atomic check name.
                    print(f"Detected trifluoro aromatic building block: {mol_smiles}")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Check for fluorination reactions
            rxn_smiles = node["metadata"]["mapped_reaction_smiles"]

            # Check for aromatic fluorination reactions
            if checker.check_reaction("Aromatic fluorination", rxn_smiles):
                has_fluorinated_aromatics = True
                findings_json["atomic_checks"]["named_reactions"].append("Aromatic fluorination")
                # Note: "Aromatic fluorination reaction" is a conceptual target for structural constraint, not an atomic check name.
                print(f"Detected aromatic fluorination reaction: {rxn_smiles}")

        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction.
            # Depth remains same when going from reaction to chemical.
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'mol' (chemical)
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Add structural constraint if any fluorinated aromatic building block or reaction was found
    if has_fluorinated_aromatics:
        # This structural constraint is met if any of the atomic checks were true.
        # The original JSON defines a 'count' constraint with operator ">= 1".
        # Since 'has_fluorinated_aromatics' becomes true if any of the conditions are met,
        # we can add the structural constraint if this flag is true.
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": [
                    "fluorinated aromatic halide building block",
                    "aromatic trifluoromethyl building block",
                    "Aromatic fluorination reaction"
                ],
                "operator": ">=",
                "value": 1
            }
        })

    return has_fluorinated_aromatics, findings_json
