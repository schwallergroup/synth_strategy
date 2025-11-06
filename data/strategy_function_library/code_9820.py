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
    This function detects if the synthesis involves a trifluoromethyl-containing substrate
    that remains unchanged throughout the synthesis.
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

    # Find the target molecule (root of the synthesis tree)
    target_mol_smiles = route["smiles"]
    print(f"Target molecule: {target_mol_smiles}")

    # Check if target contains CF3
    if not checker.check_fg("Trifluoro group", target_mol_smiles):
        print("Target molecule does not contain a trifluoromethyl group")
        return False, findings_json # Return early if target doesn't have CF3
    else:
        findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
        # Add the structural constraint for target containing CF3
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Trifluoro group"
                ],
                "comment": "The final product must contain a Trifluoro group."
            }
        })

    # Track CF3 groups through the synthesis
    cf3_preserved = True

    def dfs_traverse(node, depth=0):
        nonlocal cf3_preserved, findings_json

        if node["type"] == "reaction":
            try:
                # Get reaction SMILES with atom mapping
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has CF3
                product_has_cf3 = checker.check_fg("Trifluoro group", product)
                if product_has_cf3:
                    if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

                # Check if any reactant has CF3
                reactants_with_cf3 = [
                    r for r in reactants if checker.check_fg("Trifluoro group", r)
                ]
                if reactants_with_cf3:
                    if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

                # If product has CF3 but no reactant does, CF3 was created in this step
                if product_has_cf3 and not reactants_with_cf3:
                    print(f"CF3 group was created at depth {depth}, not preserved from substrate")
                    cf3_preserved = False
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "formation of Trifluoro group",
                            "comment": "A Trifluoro group cannot be formed during any reaction step."
                        }
                    })

                # If a reactant has CF3 but product doesn't, CF3 was modified
                if reactants_with_cf3 and not product_has_cf3:
                    print(f"CF3 group was modified/removed at depth {depth}")
                    cf3_preserved = False
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "destruction of Trifluoro group",
                            "comment": "A Trifluoro group cannot be removed or modified during any reaction step."
                        }
                    })

                print(
                    f"Reaction at depth {depth}: CF3 in product: {product_has_cf3}, CF3 in reactants: {len(reactants_with_cf3) > 0}"
                )

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

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

    # The strategy is valid if the target has CF3 and it was preserved throughout synthesis
    result = checker.check_fg("Trifluoro group", target_mol_smiles) and cf3_preserved
    print(f"Trifluoromethyl containing substrate strategy: {result}")
    return result, findings_json
