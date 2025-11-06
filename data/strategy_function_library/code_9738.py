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
    Detects a strategy where a halogen (bromine) is retained throughout the synthesis
    while another halogen (chlorine) is introduced in a late stage.
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

    # Initialize tracking variables
    depths_with_bromine = set()
    max_depth = 0
    has_late_chlorination = False

    def dfs_traverse(node, depth=0):
        nonlocal depths_with_bromine, max_depth, has_late_chlorination, findings_json

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            # Check for bromine in product
            if checker.check_fg("Aromatic halide", product_str) and "Br" in product_str:
                depths_with_bromine.add(depth)
                if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

            # Check for late-stage chlorination (depth 0 or 1)
            if depth <= 1:
                if ("Cl" in product_str and "Cl" not in reactants_str):
                    # Check for specific halide types
                    halide_found = False
                    for fg_name in ["Primary halide", "Secondary halide", "Tertiary halide", "Aromatic halide"]:
                        if checker.check_fg(fg_name, product_str):
                            if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                            halide_found = True

                    if halide_found:
                        has_late_chlorination = True
                        # Assuming 'chlorination' reaction name is detected here conceptually
                        if "chlorination" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("chlorination")

        elif node["type"] == "mol" and not node.get("in_stock", False):
            # Check for bromine in intermediate molecules
            if checker.check_fg("Aromatic halide", node["smiles"]) and "Br" in node["smiles"]:
                depths_with_bromine.add(depth)
                if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node['type'] is 'mol' or other types
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if bromine is present at all depths
    all_depths_have_bromine = len(depths_with_bromine) > 0 and all(
        d in depths_with_bromine for d in range(max_depth + 1)
    )

    # Determine final result
    result = all_depths_have_bromine and has_late_chlorination

    # Record structural constraints if conditions are met
    if all_depths_have_bromine:
        findings_json["structural_constraints"].append({
            "type": "persistence",
            "details": {
                "target": "Aromatic Bromide",
                "scope": "all_stages"
            }
        })
    if has_late_chlorination:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "chlorination",
                "position": "late_stage"
            }
        })

    # Return True if bromine is retained throughout and chlorine is introduced late
    return result, findings_json
