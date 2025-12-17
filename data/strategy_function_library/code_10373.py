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
    Detects if the synthesis follows a linear strategy (as opposed to convergent).
    Linear synthesis typically has one main reactant and one or more reagents in each step.
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

    is_linear = True

    # Track branching in the synthesis tree
    branch_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, branch_count, findings_json

        # Track branching - a convergent route is defined by a molecule being formed from multiple
        # distinct reaction steps. We check if any molecule node has more than one reaction as a child.
        if node["type"] == "mol" and not node.get("in_stock", False):
            reaction_children = [
                child for child in node.get("children", []) if child["type"] == "reaction"
            ]
            if len(reaction_children) > 1:
                branch_count += 1
                is_linear = False
                # Add the structural constraint finding
                findings_json["structural_constraints"].append({
                    "type": "count",
                    "details": {
                        "target": "convergent_steps",
                        "operator": "==",
                        "value": 0,
                        "definition": "A convergent step is defined as a molecule in the synthesis tree that is the product of more than one distinct reaction."
                    }
                })

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or 'mol' type
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return is_linear, findings_json
