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
    This function detects a synthetic strategy involving mid-stage cyanation of an aryl halide.
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

    cyanation_found = False
    max_depth = 0

    # First pass to determine max depth
    def calculate_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            calculate_max_depth(child, depth + 1)

    calculate_max_depth(route)
    print(f"Maximum depth of synthesis route: {max_depth}")

    def dfs_traverse(node, depth=0):
        nonlocal cyanation_found, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for aryl halide in reactants
                has_aryl_halide = False

                for reactant in reactants_smiles:
                    # Check for aryl halide
                    if checker.check_fg("Aromatic halide", reactant):
                        has_aryl_halide = True
                        if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                        print(f"Found aromatic halide in reactant: {reactant}")

                # Check for nitrile in product
                if checker.check_fg("Nitrile", product_smiles):
                    if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                    print(f"Found nitrile in product: {product_smiles}")

                    if has_aryl_halide:
                        # Add co-occurrence constraint if both are found in the same reaction step
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "Aromatic halide",
                                    "Nitrile"
                                ],
                                "scope": "single_reaction_step"
                            }
                        })

                        # Define mid-stage as between 30% and 70% of the maximum depth
                        mid_stage_min = max(1, int(max_depth * 0.3))
                        mid_stage_max = int(max_depth * 0.7)

                        if mid_stage_min <= depth <= mid_stage_max:
                            print(f"Detected mid-stage cyanation at depth {depth}/{max_depth}")
                            cyanation_found = True
                            if "cyanation_of_aryl_halide" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("cyanation_of_aryl_halide")
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": "cyanation_of_aryl_halide",
                                    "position": "mid_stage"
                                }
                            })
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # chemical node
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    return cyanation_found, findings_json
