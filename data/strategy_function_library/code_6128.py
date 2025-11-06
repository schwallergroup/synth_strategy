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
    This function detects if a ring formation occurs in the second half of the synthesis.
    Late stage is defined as occurring at depth <= max_depth/2.
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

    max_depth = 0
    ring_formation_depths = []
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, ring_formation_depths, findings_json
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Count rings in reactants and product
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product) if product else None

                    if product_mol and all(reactant_mols):
                        # Get accurate ring counts
                        reactant_ring_count = sum(
                            mol.GetRingInfo().NumRings() for mol in reactant_mols
                        )
                        product_ring_count = product_mol.GetRingInfo().NumRings()

                        print(
                            f"Depth {depth}: Reactant rings: {reactant_ring_count}, Product rings: {product_ring_count}"
                        )

                        # Check if this is a ring formation reaction
                        is_ring_formation = False

                        # Method 1: Check by ring count
                        if product_ring_count > reactant_ring_count:
                            is_ring_formation = True
                            print(f"Detected ring formation by count at depth {depth}")

                        if is_ring_formation:
                            ring_formation_depths.append(depth)
                            # Add to findings_json for atomic check
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                except Exception as e:
                    print(f"Error analyzing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same

            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    print(f"Maximum depth: {max_depth}")
    print(f"Ring formations detected at depths: {ring_formation_depths}")

    # Check if any ring formation occurs in the second half of synthesis
    if ring_formation_depths:
        for depth in ring_formation_depths:
            if depth <= max_depth / 2:  # Late stage (remember depth 0 is the final product)
                print(f"Late stage ring formation detected at depth {depth}")
                result = True
                # Add to findings_json for structural constraint
                if {"type": "positional", "details": {"target": "ring_formation", "position": "late_stage"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ring_formation", "position": "late_stage"}})
                break # Only need to find one to set result to True

    if not result:
        print("No late stage ring formation detected")

    return result, findings_json