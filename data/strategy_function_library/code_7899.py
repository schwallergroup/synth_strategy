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
    Detects if a synthesis route involves a borylation reaction whose product is consumed in a subsequent synthetic step along the same linear path.
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
    result = False

    # Track borylation reactions and their products
    borylation_reactions = []

    def dfs_traverse(node, depth=0, path=None):
        nonlocal result, findings_json
        if path is None:
            path = []

        current_path = path + [node]

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            product = rsmi.split(">")[-1]

            # Check for borylation reactions
            is_borylation = False

            # Check for specific borylation reaction types
            borylation_reaction_names = [
                "Preparation of boronic acids",
                "Preparation of boronic acids without boronic ether",
                "Preparation of boronic acids from trifluoroborates",
                "Preparation of boronic ethers"
            ]

            for r_name in borylation_reaction_names:
                if checker.check_reaction(r_name, rsmi):
                    is_borylation = True
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)

            # If borylation detected, store the reaction and product for later checking
            if is_borylation:
                borylation_reactions.append((product, depth, current_path))

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is chemical, increase depth
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            if dfs_traverse(child, next_depth, current_path):
                result = True # Set result to True if any path returns True
                return True

        return False

    # Start DFS traversal
    dfs_traverse(route)

    # After traversal, check if any borylation product is used in a subsequent reaction
    final_check_passed = False
    for borylation_product, depth, path in borylation_reactions:
        # Check if this product is used in any subsequent reaction
        for i, node in enumerate(path):
            if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                # Check reactions that occur AFTER this borylation (i.e., have a smaller depth)
                if i < depth:
                    # Check if the borylated product is used in this reaction
                    for reactant in reactants:
                        if reactant == borylation_product:
                            final_check_passed = True
                            result = True # Set overall result to True
                            # Add the structural constraint if it's not already there
                            structural_constraint_obj = {
                                "type": "sequence",
                                "details": {
                                    "producer_reactions": [
                                        "Preparation of boronic acids",
                                        "Preparation of boronic acids without boronic ether",
                                        "Preparation of boronic acids from trifluoroborates",
                                        "Preparation of boronic ethers"
                                    ],
                                    "consumer_reaction": "any",
                                    "relationship": "product_of_producer_is_reactant_of_consumer"
                                }
                            }
                            if structural_constraint_obj not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(structural_constraint_obj)
                            break # Found usage, no need to check other reactants for this reaction
                    if final_check_passed:
                        break # Found usage in a subsequent reaction, no need to check other reactions in this path
            if final_check_passed:
                break # Found usage in a subsequent reaction, no need to check other nodes in this path
        if final_check_passed:
            break # Found usage for one borylation product, no need to check other borylation products

    return result, findings_json