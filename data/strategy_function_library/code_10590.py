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
    This function detects a synthetic strategy where a quinoline scaffold
    is preserved throughout the synthesis.
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

    scaffold_preserved = True
    reaction_count = 0
    result = False

    # First check if the target molecule contains the quinoline scaffold
    if route["type"] == "mol" and route["smiles"]:
        target_mol = route["smiles"]
        has_quinoline = checker.check_ring("quinoline", target_mol)

        # If target doesn't have the scaffold, strategy doesn't apply
        if not has_quinoline:
            result = False
            return result, findings_json
        else:
            findings_json["atomic_checks"]["ring_systems"].append("quinoline")

    def dfs_traverse(node):
        nonlocal scaffold_preserved, reaction_count, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            product = rsmi.split(">")[-1]

            reaction_count += 1

            # Check if the quinoline scaffold is present in the product
            has_quinoline = checker.check_ring("quinoline", product)

            if not has_quinoline:
                scaffold_preserved = False
            else:
                # If quinoline is present in product, record it as detected
                if "quinoline" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("quinoline")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Only consider routes with at least 2 reactions
    result = scaffold_preserved and reaction_count >= 2

    if scaffold_preserved:
        # This corresponds to the negation constraint: "The quinoline scaffold must be present in the final product and all intermediate products."
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "quinoline_absence_in_product",
                "description": "The quinoline scaffold must be present in the final product and all intermediate products."
            }
        })
    
    if reaction_count >= 2:
        # This corresponds to the count constraint
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json
