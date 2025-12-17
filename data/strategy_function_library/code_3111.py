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
    Detects a convergent synthesis strategy that maintains a CF3 group
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

    found_convergent_step = False
    result = True # Initialize the overall result

    # Check if target molecule has CF3
    if route["type"] == "mol":
        if not checker.check_fg("Trifluoro group", route["smiles"]):
            result = False
        else:
            findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "Trifluoro group",
                    "position": "final_product"
                }
            })

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent_step, result, findings_json

        if not result: # If overall result is already False, stop early
            return False

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for convergent synthesis (multiple fragments coming together)
            if len(reactants_smiles) >= 2:
                # Check if this is a significant fragment coupling
                significant_fragments = 0
                for r in reactants_smiles:
                    mol = Chem.MolFromSmiles(r)
                    if (
                        mol and mol.GetNumHeavyAtoms() > 5
                    ):  # Consider fragments with >5 atoms significant
                        significant_fragments += 1

                if significant_fragments >= 2:
                    found_convergent_step = True
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "convergent_reaction_step",
                            "operator": ">=",
                            "value": 1,
                            "condition": "A reaction where at least two reactants have more than 5 heavy atoms."
                        }
                    })

            # Check if CF3 is maintained in this reaction
            reactant_has_cf3 = any(checker.check_fg("Trifluoro group", r) for r in reactants_smiles)
            product_has_cf3 = checker.check_fg("Trifluoro group", product_smiles)

            if reactant_has_cf3:
                # If any reactant has CF3, record it as an atomic check
                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

            # If reactant had CF3 but product doesn't, CF3 wasn't maintained
            if reactant_has_cf3 and not product_has_cf3:
                result = False # Set overall result to False
                findings_json["structural_constraints"].append({
                    "type": "negation",
                    "details": {
                        "target": "destruction of Trifluoro group"
                    }
                })
                return False # Stop this branch of traversal

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same

            if not dfs_traverse(child, new_depth):
                return False

        return True

    # Start traversal and check if CF3 is maintained throughout
    cf3_maintained = dfs_traverse(route)

    # The final result depends on both conditions
    final_result = found_convergent_step and cf3_maintained and result

    return final_result, findings_json