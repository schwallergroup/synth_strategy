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
    This function detects a strategy involving halogen-containing aromatic compounds
    maintained throughout the synthesis.
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
    has_halogen_aromatic_in_final = False
    halogen_aromatic_maintained = True
    reaction_count = 0

    # Check if the final product (root node) has an aromatic halide
    if route["type"] == "mol":
        if checker.check_fg("Aromatic halide", route["smiles"]):
            has_halogen_aromatic_in_final = True
            findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "Aromatic halide",
                    "position": "final_product"
                }
            })
        print(f"Final product: {route['smiles']}")
        print(f"Final product has halogen-aromatic: {has_halogen_aromatic_in_final}")

    def dfs_traverse(node, depth=0):
        nonlocal halogen_aromatic_maintained, reaction_count, findings_json

        # Process reaction nodes
        if node["type"] == "reaction":
            reaction_count += 1
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # In retrosynthesis, product is the starting material and reactants are the next steps
                product_has_halogen = checker.check_fg("Aromatic halide", product)
                reactants_have_halogen = any(
                    checker.check_fg("Aromatic halide", reactant) for reactant in reactants
                )

                print(f"Depth {depth}, Reaction: {rsmi}")
                print(f"  Product has halogen-aromatic: {product_has_halogen}")
                print(f"  Reactants have halogen-aromatic: {reactants_have_halogen}")

                # If product has halogen-aromatic but none of the reactants do,
                # then the halogen-aromatic wasn't maintained (it was introduced in this step)
                if product_has_halogen and not reactants_have_halogen:
                    halogen_aromatic_maintained = False
                    # This condition means the negation constraint is NOT met for this step
                    # We only add to findings_json if the constraint IS met.

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is 'chemical' (mol), depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if:
    # 1. Final product has halogen-aromatic
    # 2. Halogen-aromatic is maintained throughout synthesis
    # 3. There's at least one reaction
    strategy_present = (
        has_halogen_aromatic_in_final and halogen_aromatic_maintained and reaction_count > 0
    )

    if halogen_aromatic_maintained:
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "introduction_of_Aromatic_halide"
            }
        })
    
    if reaction_count > 0:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "any_reaction",
                "operator": ">",
                "value": 0
            }
        })

    print(
        f"Strategy assessment: has_halogen_in_final={has_halogen_aromatic_in_final}, maintained={halogen_aromatic_maintained}, reaction_count={reaction_count}"
    )

    return strategy_present, findings_json
