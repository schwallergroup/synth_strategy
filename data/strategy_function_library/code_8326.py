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
    This function detects oxazole ring formation in the early stages of synthesis.
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

    early_stage_oxazole_found = False
    max_depth = 0

    # First pass to determine the maximum depth
    def get_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            # New logic for get_max_depth to match dfs_traverse
            if node["type"] == "reaction":
                get_max_depth(child, depth)
            else:
                get_max_depth(child, depth + 1)

    get_max_depth(route)
    print(f"Maximum depth of synthesis route: {max_depth}")

    # Second pass to find oxazole formation in early stage
    def dfs_traverse(node, depth, max_depth):
        nonlocal early_stage_oxazole_found, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains oxazole
                product_has_oxazole = checker.check_ring("oxazole", product_smiles)
                if product_has_oxazole:
                    if "oxazole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("oxazole")

                # Check if any reactant contains oxazole
                reactants_have_oxazole = any(
                    checker.check_ring("oxazole", reactant) for reactant in reactants_smiles
                )

                # Check if this is an oxazole formation reaction
                if product_has_oxazole and not reactants_have_oxazole:
                    print(f"Oxazole formation detected at depth {depth}")
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                    # Check if this is early stage (depth > max_depth/2)
                    if depth > (max_depth / 2):
                        print(
                            f"Early stage oxazole formation confirmed at depth {depth} (max depth: {max_depth})"
                        )
                        early_stage_oxazole_found = True
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "oxazole_ring_formation",
                                "position": "early_stage"
                            }
                        })

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # Modified recursive call for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth, max_depth)
            else:
                dfs_traverse(child, depth + 1, max_depth)

    # Start traversal
    dfs_traverse(route, 0, max_depth)
    print(f"Early stage oxazole formation found: {early_stage_oxazole_found}")
    return early_stage_oxazole_found, findings_json
