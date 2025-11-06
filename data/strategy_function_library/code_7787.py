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


INDOLE_LIKE_SCAFFOLDS = ['indole', 'quinoline', 'isoquinoline']

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where a specific indole-like core scaffold (from the INDOLE_LIKE_SCAFFOLDS list, including indole, quinoline, and isoquinoline) is maintained. The strategy is flagged if the final product contains one of the scaffolds and the scaffold is also present in both a reactant and the product for at least three reactions in the synthesis.
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

    # Track if the final product contains a target scaffold
    final_product_has_indole = False

    # Track reactions where the scaffold core is maintained
    reactions_maintaining_indole = 0

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_indole, reactions_maintaining_indole, findings_json

        # Check if current node is a molecule
        if node["type"] == "mol" and depth == 0:
            # This is the final product (depth 0)
            mol_smiles = node["smiles"]
            # Check for a target scaffold structure
            for scaffold in INDOLE_LIKE_SCAFFOLDS:
                if checker.check_ring(scaffold, mol_smiles):
                    final_product_has_indole = True
                    if scaffold not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(scaffold)
            

        # Check if current node is a reaction
        elif node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains a target scaffold structure
            product_has_indole = False
            for scaffold in INDOLE_LIKE_SCAFFOLDS:
                if checker.check_ring(scaffold, product):
                    product_has_indole = True
                    if scaffold not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(scaffold)

            # Check if any reactant contains a target scaffold structure
            reactants_have_indole = False
            for reactant in reactants:
                for scaffold in INDOLE_LIKE_SCAFFOLDS:
                    if checker.check_ring(scaffold, reactant):
                        reactants_have_indole = True
                        if scaffold not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(scaffold)

            # If both reactants and product have the scaffold, the core is maintained
            if product_has_indole and reactants_have_indole:
                reactions_maintaining_indole += 1

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when going from chemical to reaction.
            # Depth remains the same when going from reaction to chemical.
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'mol' (chemical)
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    result = final_product_has_indole and reactions_maintaining_indole >= 3

    # Add structural constraints if conditions are met
    if final_product_has_indole:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": [
                    "indole",
                    "quinoline",
                    "isoquinoline"
                ],
                "position": "last_stage"
            }
        })
    
    if reactions_maintaining_indole >= 3:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "indole_like_scaffold_maintenance",
                "operator": ">=",
                "value": 3
            }
        })

    # Return True if:
    # 1. The final product has the scaffold
    # 2. At least 3 reactions maintain the scaffold core
    return result, findings_json
