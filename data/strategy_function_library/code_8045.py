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
    Detects if a cyclopropane motif is present in the final product and maintained
    throughout the synthesis.
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

    # Check if final product has cyclopropane
    final_product_has_cyclopropane = False
    result = False

    if route["type"] == "mol":
        final_product_has_cyclopropane = checker.check_ring("cyclopropane", route["smiles"])
        print(f"Final product has cyclopropane: {final_product_has_cyclopropane}")

        if final_product_has_cyclopropane:
            findings_json["atomic_checks"]["ring_systems"].append("cyclopropane")
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "cyclopropane",
                    "position": "final_product"
                }
            })
        else:
            return False, findings_json
    else:
        print("Error: Root node is not a molecule")
        return False, findings_json

    # Track if cyclopropane is maintained throughout synthesis
    cyclopropane_broken = False

    def dfs_traverse(node, depth=0):
        nonlocal cyclopropane_broken, findings_json

        if cyclopropane_broken:
            return  # Stop traversal if we already found a break

        if node["type"] == "reaction":
            # Check if reaction maintains cyclopropane
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                product = rsmi.split(">")[-1]
                reactants = rsmi.split(">")[0].split(".")

                product_has_cyclopropane = checker.check_ring("cyclopropane", product)

                # Check which reactants have cyclopropane
                reactants_with_cyclopropane = [
                    r for r in reactants if checker.check_ring("cyclopropane", r)
                ]
                reactant_has_cyclopropane = len(reactants_with_cyclopropane) > 0

                print(
                    f"Reaction at depth {depth}: product has cyclopropane: {product_has_cyclopropane}, reactants with cyclopropane: {len(reactants_with_cyclopropane)}"
                )

                # In retrosynthesis:
                # If product has cyclopropane but no reactant does, cyclopropane is formed in this step (forward)
                # This is fine - we're looking for maintenance, not formation
                if product_has_cyclopropane and not reactant_has_cyclopropane:
                    print(f"Cyclopropane formed in reaction at depth {depth} (forward direction)")
                    # This implies a 'ring_formation' reaction for cyclopropane
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                # If product doesn't have cyclopropane but a reactant does, it's broken in this step (forward)
                elif not product_has_cyclopropane and reactant_has_cyclopropane:
                    cyclopropane_broken = True
                    print(f"Cyclopropane broken in reaction at depth {depth} (forward direction)")
                    # This implies a 'ring_destruction' reaction for cyclopropane
                    if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # If final product has cyclopropane and it's never broken, return True
    result = final_product_has_cyclopropane and not cyclopropane_broken
    print(f"Cyclopropane maintained: {result}")

    if not cyclopropane_broken:
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "ring_destruction",
                "context": {
                    "ring_system": "cyclopropane"
                }
            }
        })

    return result, findings_json
