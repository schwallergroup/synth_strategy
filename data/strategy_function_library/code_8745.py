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
    This function detects if the synthesis uses a late-stage cyanation strategy,
    where a nitrile group is introduced in the final or penultimate step.
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

    final_product_has_nitrile = False
    cyanation_at_depth = None
    final_product_smiles = None

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_nitrile, cyanation_at_depth, final_product_smiles, findings_json

        if node["type"] == "mol" and depth == 0:  # Final product
            final_product_smiles = node["smiles"]
            if checker.check_fg("Nitrile", final_product_smiles):
                final_product_has_nitrile = True
                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                # Record structural constraint for Nitrile in final product
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "Nitrile",
                        "position": "final_product"
                    }
                })

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # A robust check for cyanation is a net increase in nitrile groups.
                product_nitrile_count = len(checker.get_fg_atom_indices("Nitrile", product))
                reactants_nitrile_count = sum(
                    len(checker.get_fg_atom_indices("Nitrile", r)) for r in reactants
                )

                if product_nitrile_count > reactants_nitrile_count:
                    # This implies a cyanation reaction occurred
                    findings_json["atomic_checks"]["named_reactions"].append("cyanation")
                    if cyanation_at_depth is None or depth < cyanation_at_depth:
                        cyanation_at_depth = depth

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., chemical), increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Late stage is defined as depth 0 or 1
    is_late_stage = cyanation_at_depth is not None and cyanation_at_depth <= 1

    # The function should return True only if a cyanation was *actually detected* at a late stage
    # and the final product contains a nitrile.
    result = final_product_has_nitrile and is_late_stage

    if is_late_stage:
        # Record structural constraint for cyanation at late stage
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "cyanation",
                "position_constraint": {
                    "operator": "<=",
                    "value": 1,
                    "unit": "depth"
                }
            }
        })

    return result, findings_json
