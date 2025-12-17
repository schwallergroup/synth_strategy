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
    Detects if an aromatic halide pattern is maintained throughout the synthesis.
    This means that once an aromatic halide group (Ar-F, Ar-Cl, Ar-Br, Ar-I) is present in an intermediate,
    it is preserved in all subsequent intermediates leading to the final product.
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

    result = True # Initialize the main boolean result

    # Track molecules with chloroaromatic groups
    molecules_with_chloroaromatic = {}

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Skip starting materials (in_stock)
            if node.get("in_stock", False):
                return

            # Check for chloroaromatic pattern using the checker
            has_chloroaromatic = checker.check_fg("Aromatic halide", mol_smiles)
            if has_chloroaromatic:
                # Add to atomic checks if found
                if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

            # Store result with depth information
            molecules_with_chloroaromatic[mol_smiles] = {
                "has_chloroaromatic": has_chloroaromatic,
                "depth": depth,
            }

        # Traverse children (going deeper in retrosynthesis)
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'mol' (chemical)
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # If no non-starting materials were found, return False
    if not molecules_with_chloroaromatic:
        result = False
        return result, findings_json

    # Check if the final product (depth 0) has a chloroaromatic group
    final_product = None
    for mol_smiles, info in molecules_with_chloroaromatic.items():
        if info["depth"] == 0:
            final_product = mol_smiles
            break

    if final_product is None:
        result = False
        return result, findings_json

    if not molecules_with_chloroaromatic[final_product]["has_chloroaromatic"]:
        result = False
        return result, findings_json
    else:
        # If final product has Aromatic halide, add positional constraint
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Aromatic halide",
                "position": "final_product"
            }
        })

    # Check if chloroaromatic pattern is maintained throughout the synthesis
    # This means if a molecule has a chloroaromatic group, all molecules with lower depth
    # (closer to the final product) should also have it
    for mol_smiles, info in molecules_with_chloroaromatic.items():
        if info["has_chloroaromatic"]:
            # Check all molecules with lower depth
            for other_smiles, other_info in molecules_with_chloroaromatic.items():
                if other_info["depth"] < info["depth"]:
                    if not other_info["has_chloroaromatic"]:
                        result = False
                        # If a break in the pattern is found, add the negation constraint
                        findings_json["structural_constraints"].append({
                            "type": "negation",
                            "details": {
                                "target": "sequence",
                                "description": "A molecule containing an 'Aromatic halide' is followed in a later synthesis stage (closer to the product) by a molecule that does not contain it."
                            }
                        })
                        return result, findings_json

    return result, findings_json