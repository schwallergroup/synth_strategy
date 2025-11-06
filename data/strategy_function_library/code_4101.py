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
    Detects a synthesis strategy centered around a piperazine scaffold.

    A piperazine scaffold strategy involves:
    1. Formation or use of a piperazine ring
    2. The piperazine ring is central to the synthesis (appears in multiple steps)
    3. The final product contains a piperazine ring
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

    # Track piperazine occurrences and their depths
    piperazine_nodes = []
    max_depth = 0
    final_product_has_piperazine = False
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, final_product_has_piperazine, findings_json

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if smiles:
                # Check if molecule contains piperazine
                if checker.check_ring("piperazine", smiles):
                    piperazine_nodes.append((depth, "mol", smiles))
                    if "piperazine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("piperazine")

                    # If this is the final product (depth 0), mark it
                    if depth == 0:
                        final_product_has_piperazine = True

        elif node["type"] == "reaction":
            # Check if this reaction forms or modifies a piperazine ring
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has piperazine but reactants don't (piperazine formation)
                product_has_piperazine = checker.check_ring("piperazine", product)
                reactants_have_piperazine = any(
                    checker.check_ring("piperazine", r) for r in reactants
                )

                if product_has_piperazine:
                    if "piperazine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("piperazine")
                    if not reactants_have_piperazine:
                        piperazine_nodes.append((depth, "formation", rsmi))
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                    else:
                        piperazine_nodes.append((depth, "modification", rsmi))

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Analyze results to determine if this is a piperazine scaffold strategy
    if not piperazine_nodes:
        return False, findings_json

    # Check if final product has piperazine
    if not final_product_has_piperazine:
        return False, findings_json
    else:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "piperazine",
                "position": "last_stage"
            }
        })

    # Check if piperazine appears in multiple steps
    unique_depths = len(set(node[0] for node in piperazine_nodes))
    if unique_depths < 2:
        return False, findings_json
    else:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "piperazine_presence_steps",
                "operator": ">=",
                "value": 2
            }
        })

    # Check if piperazine is formed early
    formation_nodes = [node for node in piperazine_nodes if node[1] == "formation"]
    if formation_nodes:
        earliest_formation = min(node[0] for node in formation_nodes)
        # A late-stage formation is defined as occurring in the last 30% of steps.
        # A low depth value corresponds to a late stage.
        if earliest_formation < (max_depth * 0.3):
            # This condition means it was formed late, which is a NEGATION of the desired strategy
            # So, if this is true, the strategy fails.
            findings_json["structural_constraints"].append({
                "type": "negation",
                "details": {
                    "type": "positional",
                    "details": {
                        "target": "piperazine_formation",
                        "position": "late_stage"
                    }
                }
            })
            return False, findings_json

    result = True
    return result, findings_json
