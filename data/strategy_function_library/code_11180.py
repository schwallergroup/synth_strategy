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


HETEROCYCLES_OF_INTEREST = [
    "furan",
    "pyrrole",
    "thiophene",
    "pyrazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "triazole",
    "tetrazole",
    "pyridine",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "indole",
    "benzofuran",
    "benzothiophene",
    "benzimidazole",
    "benzoxazole",
    "benzothiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis uses a trifluoromethyl-substituted heterocycle as a key building block
    that persists through to the final product.
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

    # Track molecules with CF3-heterocycle at each depth
    cf3_heterocycles = {}

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_cf3 = checker.check_fg("Trifluoro group", mol_smiles)
            if has_cf3:
                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

            # Check for heterocycles
            heterocycle_found = None
            for het_type in HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(het_type, mol_smiles):
                    heterocycle_found = het_type
                    if het_type not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(het_type)
                    break

            # If molecule has both CF3 and heterocycle
            if has_cf3 and heterocycle_found:
                # Record the co-occurrence constraint if not already recorded
                co_occurrence_constraint = {
                    "type": "co-occurrence",
                    "details": {
                        "targets": [
                            "Trifluoro group",
                            "any_specified_heterocycle"
                        ],
                        "scope": "per_molecule",
                        "description": "A key building block is defined as a single molecule containing both a Trifluoro group and one of the heterocycles listed in atomic_checks."
                    }
                }
                if co_occurrence_constraint not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(co_occurrence_constraint)

                # Store this molecule at its depth
                if depth not in cf3_heterocycles:
                    cf3_heterocycles[depth] = []
                cf3_heterocycles[depth].append(
                    (mol_smiles, heterocycle_found, node.get("in_stock", False))
                )

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is 'mol', depth increases
            next_depth = depth + 1

        # Traverse children (moving backward in synthesis)
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the final product (root node)
    dfs_traverse(route)

    # Check if CF3-heterocycle persists from starting material to final product
    valid_strategy = False

    if cf3_heterocycles:
        # Check if we have CF3-heterocycle at depth 0 (final product)
        has_cf3_in_final_product = False
        if 0 in cf3_heterocycles:
            has_cf3_in_final_product = True

        # Check if we have CF3-heterocycle in any starting material
        starting_materials_with_cf3 = []
        for depth in cf3_heterocycles:
            for smiles, het_type, in_stock in cf3_heterocycles[depth]:
                if in_stock:
                    starting_materials_with_cf3.append((depth, smiles, het_type))
        has_cf3_in_starting_material = len(starting_materials_with_cf3) > 0

        # Check if we have at least one intermediate with CF3-heterocycle
        has_intermediate = any(depth > 0 for depth in cf3_heterocycles.keys() if depth != 0)

        # If we have CF3-heterocycle in final product, at least one intermediate,
        # and at least one starting material, then the strategy is valid
        if has_cf3_in_final_product and has_intermediate and has_cf3_in_starting_material:
            valid_strategy = True
            # Record the persistence constraint
            persistence_constraint = {
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "key_building_block_in_starting_material",
                        "key_building_block_in_intermediate",
                        "key_building_block_in_final_product"
                    ],
                    "scope": "per_route",
                    "description": "The key building block must be present in a starting material, at least one intermediate, and the final product, indicating its persistence throughout the synthesis."
                }
            }
            if persistence_constraint not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append(persistence_constraint)

    return valid_strategy, findings_json