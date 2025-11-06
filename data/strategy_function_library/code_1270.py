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
    This function detects the formation of a quinazoline core during the synthesis.
    """
    from rdkit import Chem

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    quinazoline_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal quinazoline_formation_detected, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains a quinazoline structure
                try:
                    product_has_quinazoline = checker.check_ring("quinazoline", product_smiles)

                    if product_has_quinazoline:
                        findings_json["atomic_checks"]["ring_systems"].append("quinazoline")
                        print(f"Product contains quinazoline at depth {depth}: {product_smiles}")

                        # Check if any reactant has quinazoline
                        reactant_has_quinazoline = False
                        for reactant_smiles in reactants_smiles:
                            if checker.check_ring("quinazoline", reactant_smiles):
                                reactant_has_quinazoline = True
                                print(f"Reactant already contains quinazoline: {reactant_smiles}")
                                break

                        if not reactant_has_quinazoline:
                            print(
                                f"Quinazoline core formation detected in reaction at depth {depth}: {rsmi}"
                            )
                            quinazoline_formation_detected = True
                except Exception as e:
                    print(f"Error checking quinazoline structure: {e}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)
    print(f"Quinazoline formation detected: {quinazoline_formation_detected}")

    return quinazoline_formation_detected, findings_json
