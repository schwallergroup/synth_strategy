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
    Detects a synthetic strategy that maintains a cyclopropane ring throughout
    the synthesis without modification once it's introduced.
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

    cyclopropane_introduced = False
    cyclopropane_preserved = True

    def dfs_traverse(node, depth=0):
        nonlocal cyclopropane_introduced, cyclopropane_preserved, findings_json

        if not cyclopropane_preserved:
            return  # Early termination if we already found a violation

        if node["type"] == "mol" and "smiles" in node:
            # Check if molecule contains a cyclopropane ring
            has_cyclopropane = checker.check_ring("cyclopropane", node["smiles"])

            if has_cyclopropane:
                cyclopropane_introduced = True
                findings_json["atomic_checks"]["ring_systems"].append("cyclopropane")
                # print(f"Molecule at depth {depth} has cyclopropane: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Only check preservation if cyclopropane has been introduced
            if cyclopropane_introduced:
                try:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if cyclopropane exists in product and reactants
                    product_has_cyclopropane = checker.check_ring("cyclopropane", product)
                    reactants_with_cyclopropane = [
                        r for r in reactants if checker.check_ring("cyclopropane", r)
                    ]

                    # If any reactant has cyclopropane but product doesn't, cyclopropane was destroyed
                    if reactants_with_cyclopropane and not product_has_cyclopropane:
                        # print(f"Reaction at depth {depth} destroys cyclopropane: {rsmi}")
                        cyclopropane_preserved = False
                        findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")
                except Exception as e:
                    # print(f"Error analyzing reaction at depth {depth}: {e}")
                    pass # Suppress print for refactored code

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    result = cyclopropane_introduced and cyclopropane_preserved

    # Add structural constraints based on final flags
    if cyclopropane_introduced:
        # This corresponds to the 'count' constraint for cyclopropane
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "ring_system:cyclopropane",
                "operator": ">=",
                "value": 1
            }
        })
    if cyclopropane_introduced and cyclopropane_preserved:
        # This corresponds to the 'negation' constraint for ring_destruction
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "ring_destruction",
                "scope": "cyclopropane"
            }
        })

    # Return True if cyclopropane was introduced and preserved throughout
    return result, findings_json
