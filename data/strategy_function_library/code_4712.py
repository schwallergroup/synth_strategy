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
    Detects if the synthesis route incorporates a trifluoromethyl (CF3) group
    that is preserved in the final product.
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

    # Track if CF3 is in final product and if it was incorporated during synthesis
    cf3_in_final = False
    cf3_incorporated = False

    # Track molecules with CF3 groups by their SMILES
    cf3_containing_mols = set()
    starting_materials_with_cf3 = set()

    def dfs_traverse(node, depth=0):
        nonlocal cf3_in_final, cf3_incorporated, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if molecule contains CF3 group
            has_cf3 = checker.check_fg("Trifluoro group", mol_smiles)

            if has_cf3:
                print(f"Detected CF3 group in molecule at depth {depth}: {mol_smiles}")
                findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
                cf3_containing_mols.add(mol_smiles)

                # If this is a starting material with CF3
                if node.get("in_stock", False):
                    starting_materials_with_cf3.add(mol_smiles)
                    print(f"Starting material contains CF3: {mol_smiles}")

                # If this is the final product
                if depth == 0:
                    cf3_in_final = True
                    print("Final product contains CF3 group")

        # Determine the new depth for recursive calls
        new_depth = depth
        if node["type"] != "reaction":  # Only increase depth if not a reaction node
            new_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, new_depth) 

    # Traverse the synthesis route
    dfs_traverse(route)

    result = False
    # Check if CF3 was incorporated during synthesis and preserved to final product
    if cf3_in_final:
        # If final product has CF3 and either:
        # 1. CF3 was directly incorporated in a reaction step, or
        # 2. None of the starting materials had CF3 (meaning it must have been created during synthesis)
        if cf3_incorporated or (
            len(starting_materials_with_cf3) == 0 and len(cf3_containing_mols) > 0
        ):
            print("CF3 group was incorporated during synthesis and preserved in final product")
            result = True
            # Add structural constraints if the condition is met
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "Trifluoro group",
                    "position": "last_stage"
                }
            })
            if len(starting_materials_with_cf3) == 0:
                findings_json["structural_constraints"].append({
                    "type": "negation",
                    "details": {
                        "target": "Trifluoro group",
                        "position": "starting_material"
                    }
                })
        else:
            print(
                "CF3 group was present in starting materials but not incorporated during synthesis"
            )
    else:
        print("Final product does not contain a CF3 group")

    return result, findings_json
