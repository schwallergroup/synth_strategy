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
    This function detects a synthetic strategy involving late-stage ester hydrolysis
    to generate a carboxylic acid as the final product.
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

    found_ester_hydrolysis = False
    final_product_has_acid = False

    # First check if the final product (root node) has a carboxylic acid
    if route["type"] == "mol" and checker.check_fg("Carboxylic acid", route["smiles"]):
        final_product_has_acid = True
        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Carboxylic acid",
                "position": "final_product"
            }
        })
        print(f"Final product has carboxylic acid: {route['smiles']}")
    else:
        print(f"Final product does not have carboxylic acid or is not a molecule")
        return False, findings_json

    def dfs_traverse(node, depth=0):
        nonlocal found_ester_hydrolysis, findings_json

        # Check for reaction nodes
        if node["type"] == "reaction":
            # Depth 1 means it's the final reaction step (directly connected to final product)
            if depth == 1:
                print(f"Examining reaction at depth {depth}")

                try:
                    # Get reaction SMILES
                    rsmi = node.get("metadata", {}).get("rsmi", "")
                    if not rsmi:
                        print("No reaction SMILES found in final step")
                        return

                    # Parse reaction components
                    parts = rsmi.split(">")
                    if len(parts) < 3:
                        print(f"Invalid reaction SMILES format: {rsmi}")
                        return

                    reactants = parts[0].split(".")
                    product = parts[2]

                    print(f"Checking reaction: {rsmi}")
                    print(f"Reactants: {reactants}")
                    print(f"Product: {product}")

                    # Check for ester in reactants
                    has_ester = False
                    for r in reactants:
                        if r and checker.check_fg("Ester", r):
                            has_ester = True
                            if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Ester")
                            break

                    # Check for carboxylic acid in product
                    has_acid = checker.check_fg("Carboxylic acid", product)
                    if has_acid and "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                    print(f"Has ester in reactants: {has_ester}")
                    print(f"Has carboxylic acid in product: {has_acid}")

                    if has_ester and has_acid:
                        # Add co-occurrence constraint if both are found
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "Ester_in_reactant",
                                    "Carboxylic acid_in_product"
                                ],
                                "scope": "last_stage_reaction"
                            }
                        })

                        # Verify this is an ester hydrolysis reaction using multiple reaction checks
                        is_hydrolysis = False
                        hydrolysis_reactions = [
                            "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                            "Ester saponification (methyl deprotection)",
                            "Ester saponification (alkyl deprotection)"
                        ]
                        for r_name in hydrolysis_reactions:
                            if checker.check_reaction(r_name, rsmi):
                                is_hydrolysis = True
                                if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append(r_name)

                        print(f"Is hydrolysis reaction: {is_hydrolysis}")

                        if is_hydrolysis:
                            found_ester_hydrolysis = True
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": "ester_hydrolysis",
                                    "position": "last_stage"
                                }
                            })
                            print("Found ester hydrolysis in final reaction step")

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Continue traversing the tree
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if ester hydrolysis is found in the final reaction step
    # and the final product has a carboxylic acid
    strategy_found = found_ester_hydrolysis and final_product_has_acid
    print(f"Late-stage ester hydrolysis strategy detected: {strategy_found}")
    return strategy_found, findings_json
