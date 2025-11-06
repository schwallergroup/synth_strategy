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
    Detects a synthetic strategy where key functional groups (amide, aromatic halide)
    are maintained throughout the synthesis.
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

    # Track if amide and fluorophenyl groups are present at different stages
    final_product_has_amide = False
    final_product_has_fluorophenyl = False
    early_intermediate_has_amide = False
    early_intermediate_has_fluorophenyl = False

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_amide, final_product_has_fluorophenyl
        nonlocal early_intermediate_has_amide, early_intermediate_has_fluorophenyl
        nonlocal findings_json

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol_smiles = node["smiles"]

                current_has_amide = False
                current_has_fluorophenyl = False

                # Check for amide groups using checker functions for all amide types
                amide_types = ["Primary amide", "Secondary amide", "Tertiary amide"]
                for amide_type in amide_types:
                    if checker.check_fg(amide_type, mol_smiles):
                        current_has_amide = True
                        if amide_type not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(amide_type)

                # Check for fluorophenyl groups - ensure F is on an aromatic ring
                if checker.check_fg("Aromatic halide", mol_smiles):
                    current_has_fluorophenyl = True
                    if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

                print(
                    f"Depth {depth}: {mol_smiles} - Amide: {current_has_amide}, Fluorophenyl: {current_has_fluorophenyl}"
                )

                # Check if this is the final product (depth 0 in retrosynthesis)
                if depth == 0:
                    print(f"Final product detected at depth {depth}")
                    final_product_has_amide = current_has_amide
                    final_product_has_fluorophenyl = current_has_fluorophenyl

                # Check if this is an intermediate (depth >= 1 in retrosynthesis)
                elif depth >= 1:
                    print(f"Intermediate detected at depth {depth}")
                    if current_has_amide:
                        early_intermediate_has_amide = True
                    if current_has_fluorophenyl:
                        early_intermediate_has_fluorophenyl = True
            except Exception as e:
                print(f"Error processing molecule: {e}")

        # Continue traversing
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if both functional groups persist throughout synthesis
    has_persistent_functional_groups = (
        final_product_has_amide
        and final_product_has_fluorophenyl
        and early_intermediate_has_amide
        and early_intermediate_has_fluorophenyl
    )

    print(f"Final product has amide: {final_product_has_amide}")
    print(f"Final product has fluorophenyl: {final_product_has_fluorophenyl}")
    print(f"Early intermediate has amide: {early_intermediate_has_amide}")
    print(f"Early intermediate has fluorophenyl: {early_intermediate_has_fluorophenyl}")

    if has_persistent_functional_groups:
        print("Detected persistence of key functional groups (amide, fluorophenyl)")

    # Populate structural constraints based on the flags
    if final_product_has_amide:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": [
                    "Primary amide",
                    "Secondary amide",
                    "Tertiary amide"
                ],
                "position": "last_stage",
                "quantifier": "any"
            }
        })
    if final_product_has_fluorophenyl:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Aromatic halide",
                "position": "last_stage",
                "quantifier": "any"
            }
        })
    if early_intermediate_has_amide:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": [
                    "Primary amide",
                    "Secondary amide",
                    "Tertiary amide"
                ],
                "position": "not_last_stage",
                "quantifier": "any"
            }
        })
    if early_intermediate_has_fluorophenyl:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Aromatic halide",
                "position": "not_last_stage",
                "quantifier": "any"
            }
        })

    return has_persistent_functional_groups, findings_json
