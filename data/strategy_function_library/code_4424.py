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
    """Detects if an imidazole ring and a nitro group are concurrently present in molecules throughout a synthesis. It verifies that product molecules in each step contain both an imidazole ring and a nitro group, and that these are not lost in the corresponding reactants. This function does not validate that the nitro group is attached to the imidazole ring."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Track the presence of the functional groups
    scaffold_components_present = False
    scaffold_components_maintained = True
    reaction_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal scaffold_components_present, scaffold_components_maintained, reaction_count, findings_json

        if node["type"] == "reaction":
            reaction_count += 1
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0]
                products_smiles = rsmi.split(">")[-1]

                try:
                    product_has_imidazole = checker.check_ring("imidazole", products_smiles)
                    if product_has_imidazole:
                        if "imidazole" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("imidazole")

                    product_has_nitro = checker.check_fg("Nitro group", products_smiles)
                    if product_has_nitro:
                        if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

                    product_has_components = product_has_imidazole and product_has_nitro

                    reactant_has_imidazole = checker.check_ring("imidazole", reactants_smiles)
                    if reactant_has_imidazole:
                        if "imidazole" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("imidazole")

                    reactant_has_nitro = checker.check_fg("Nitro group", reactants_smiles)
                    if reactant_has_nitro:
                        if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

                    reactant_has_components = reactant_has_imidazole and reactant_has_nitro

                    if product_has_components:
                        scaffold_components_present = True
                        # Structural constraint: At least one reaction product in the route must contain both an imidazole ring and a Nitro group.
                        constraint_obj = {
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "imidazole",
                                    "Nitro group"
                                ],
                                "scope": "product",
                                "quantifier": "at_least_one_step",
                                "description": "At least one reaction product in the route must contain both an imidazole ring and a Nitro group."
                            }
                        }
                        if constraint_obj not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(constraint_obj)

                    if product_has_components and not reactant_has_components:
                        scaffold_components_maintained = False
                        # Structural constraint: For any reaction step, it is forbidden for the product to contain both the imidazole and nitro group if the corresponding reactant does not.
                        constraint_obj = {
                            "type": "negation",
                            "details": {
                                "description": "For any reaction step, it is forbidden for the product to contain both the imidazole and nitro group if the corresponding reactant does not."
                            }
                        }
                        if constraint_obj not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(constraint_obj)

                except Exception:
                    pass

        # Determine the new depth based on the current node's type
        new_depth = depth
        if node["type"] != "reaction":  # If current node is 'chemical' or other non-reaction type
            new_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if components are present and maintained through at least two reactions
    strategy_present = scaffold_components_present and scaffold_components_maintained and reaction_count >= 2

    # Structural constraint: The strategy must be upheld for a route with at least two reaction steps.
    if reaction_count >= 2:
        constraint_obj = {
            "type": "count",
            "details": {
                "target": "reaction",
                "operator": ">=",
                "value": 2
            }
        }
        if constraint_obj not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(constraint_obj)

    return strategy_present, findings_json