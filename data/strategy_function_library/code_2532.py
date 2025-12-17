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
    This function detects if an aromatic fluorine is retained throughout the synthesis.
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

    # Track if an aromatic fluorine is present in the final product and never removed.
    aromatic_fluorine_retained = False

    def has_aromatic_fluorine(smiles):
        """Check if a molecule has an aromatic carbon with fluorine attached"""
        # Use the specific checker for robustness and efficiency.
        if checker.check_fg("Aromatic fluorine", smiles):
            if "Aromatic fluorine" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Aromatic fluorine")
            return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_fluorine_retained, findings_json

        if node["type"] == "mol":
            # If this is the target molecule (depth 0), it must have the group
            # for the strategy to be valid.
            # This sets the initial condition.
            if depth == 0 and has_aromatic_fluorine(node["smiles"]):
                aromatic_fluorine_retained = True
                # Record the positional constraint for the final product
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "Aromatic fluorine",
                        "position": "last_stage"
                    }
                })

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            product_has_aromatic_f = has_aromatic_fluorine(product)
            reactants_have_aromatic_f = any(has_aromatic_fluorine(r) for r in reactants)

            # If the aromatic fluorine status changes (is introduced or removed),
            # then it was not retained throughout the synthesis.
            if product_has_aromatic_f != reactants_have_aromatic_f:
                aromatic_fluorine_retained = False
                if product_has_aromatic_f and not reactants_have_aromatic_f:
                    # Aromatic fluorine was introduced
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "introduction",
                            "scope": "Aromatic fluorine"
                        }
                    })
                elif not product_has_aromatic_f and reactants_have_aromatic_f:
                    # Aromatic fluorine was removed
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "removal",
                            "scope": "Aromatic fluorine"
                        }
                    })

        # Determine the new depth based on the current node's type
        new_depth = depth
        if node["type"] != "reaction":  # This means it's a 'mol' node
            new_depth = depth + 1

        # Traverse children with the calculated new depth
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # If aromatic_fluorine_retained is still True after traversal, it means
    # it was present in the final product and never introduced/removed.
    # The negation constraints are only added if the condition is violated.
    # So, if it's true, it implies the negation constraints were NOT met.
    # The positional constraint is added when the final product is checked.

    return aromatic_fluorine_retained, findings_json
