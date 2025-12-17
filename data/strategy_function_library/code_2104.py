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


BROMINATION_REACTIONS_OF_INTEREST = [
    "Aromatic bromination",
    "Bromination",
    "Wohl-Ziegler bromination benzyl primary",
    "Wohl-Ziegler bromination benzyl secondary",
    "Wohl-Ziegler bromination benzyl tertiary",
    "Wohl-Ziegler bromination allyl primary",
    "Wohl-Ziegler bromination allyl secondary",
    "Wohl-Ziegler bromination allyl tertiary",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving a late-stage Rosenmund-von Braun cyanation.
    The strategy also requires at least one preceding bromination step. If a nitro
    reduction step is present, it must occur earlier in the synthesis than the cyanation.
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

    has_cyanation = False
    has_bromination = False
    has_nitro_reduction = False
    cyanation_depth = None
    bromination_depths = []
    nitro_reduction_depth = None

    def dfs_traverse(node, depth=1):
        nonlocal has_cyanation, has_bromination, has_nitro_reduction
        nonlocal cyanation_depth, bromination_depths, nitro_reduction_depth, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for cyanation (halide to nitrile)
            if checker.check_reaction("Rosenmund-von Braun reaction", rsmi):
                has_cyanation = True
                cyanation_depth = depth
                if "Rosenmund-von Braun reaction" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Rosenmund-von Braun reaction")

            # Check for bromination reactions
            for rxn in BROMINATION_REACTIONS_OF_INTEREST:
                if checker.check_reaction(rxn, rsmi):
                    has_bromination = True
                    bromination_depths.append(depth)
                    if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)

            # Check for nitro reduction (NO2 to NH2)
            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                has_nitro_reduction = True
                nitro_reduction_depth = depth
                if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # The traversal starts at the root product, so the first reaction is at depth 1.
    # We initiate the traversal on the children of the root node.
    for child in route.get("children", []):
        dfs_traverse(child, depth=1)

    # Check if the strategy criteria are met
    strategy_present = (
        has_cyanation
        and has_bromination
        and len(bromination_depths) >= 1
        and cyanation_depth == 1  # Cyanation must be the final step
    )

    # Record structural constraints if met
    if has_cyanation and has_bromination:
        # This corresponds to the first structural constraint
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Rosenmund-von Braun reaction",
                    "bromination"
                ],
                "note": "The route must contain at least one Rosenmund-von Braun reaction and at least one reaction from the bromination list."
            }
        })

    if cyanation_depth == 1:
        # This corresponds to the second structural constraint
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Rosenmund-von Braun reaction",
                "position": "last_stage"
            }
        })

    # If nitro reduction is present, ensure it occurs earlier than cyanation
    if has_nitro_reduction:
        if nitro_reduction_depth is not None and cyanation_depth is not None and nitro_reduction_depth > cyanation_depth:
            strategy_present = strategy_present and True # Condition met
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "first": "Reduction of nitro groups to amines",
                    "second": "Rosenmund-von Braun reaction",
                    "condition": "if_present",
                    "note": "If 'Reduction of nitro groups to amines' is present, it must occur in an earlier synthetic step than the 'Rosenmund-von Braun reaction'."
                }
            })
        else:
            strategy_present = False # Condition not met

    return strategy_present, findings_json