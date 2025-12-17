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
    Detects if an oxime functional group, once present, is preserved in all subsequent reaction steps. The strategy requires the final target to contain an oxime.
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
    result = True

    # This strategy is only applicable if the final product contains an oxime.
    if not checker.check_fg("Oxime", route["smiles"]):
        result = False
    else:
        findings_json["atomic_checks"]["functional_groups"].append("Oxime")
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Oxime in final product"]}})

    def dfs_traverse(node, depth=0):
        nonlocal result, findings_json
        if not result: # Early exit if strategy already failed
            return False

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                # This parsing handles both R>>P and R>A>P formats.
                reactants_smi_list = rsmi.split(">")[0].split(".")
                product_smi = rsmi.split(">")[-1]

                reactant_has_oxime = False
                for r_smi in reactants_smi_list:
                    if checker.check_fg("Oxime", r_smi):
                        reactant_has_oxime = True
                        # No need to add to findings_json here, as it's about preservation, not initial presence
                        break

                # If an oxime is present in any reactant, it must be preserved in the product.
                if reactant_has_oxime and not checker.check_fg("Oxime", product_smi):
                    result = False  # Strategy fails: oxime was not maintained.
                    findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "destruction of Oxime"}})
                    return False
            
            except (KeyError, IndexError):
                # Fail safely if reaction data is malformed.
                result = False
                return False

        # Determine the depth for the recursive call based on the current node's type.
        if node["type"] == "reaction":
            next_depth = depth  # Depth remains the same when going from reaction to chemical
        else:
            next_depth = depth + 1  # Depth increases when going from chemical to reaction

        # Continue traversal for all children.
        for child in node.get("children", []):
            if not dfs_traverse(child, next_depth):
                return False # Propagate failure up

        return True

    if result: # Only traverse if initial check passed
        dfs_traverse(route)

    return result, findings_json