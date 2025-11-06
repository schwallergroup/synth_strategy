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


# Refactored lists of named reactions
ETHERIFICATION_REACTIONS = [
    "Williamson Ether Synthesis",
    "Williamson Ether Synthesis (intra to epoxy)",
    "Alcohol to ether",
]

ALCOHOL_HALOGENATION_REACTIONS = [
    "Alkyl iodides from alcohols",
    "Alkyl bromides from alcohols",
    "Alkyl chlorides from alcohols",
    "Appel reaction",
]

ACID_ESTER_REDUCTION_REACTIONS = [
    "Reduction of carboxylic acid to primary alcohol",
    "Reduction of ester to primary alcohol",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a specific multi-step functional group interconversion strategy. The strategy must contain the following
    sequence of reaction types in order from earliest to latest step: (1) an acid or ester reduction,
    (2) an alcohol to alkyl-halide conversion, and (3) an ether synthesis. The function validates this by checking
    for specific named reactions from predefined lists for each step: `ACID_ESTER_REDUCTION_REACTIONS`,
    `ALCOHOL_HALOGENATION_REACTIONS`, and `ETHERIFICATION_REACTIONS`.
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

    # We'll track paths that contain our target transformations
    valid_paths = []

    def dfs_traverse(node, path=None, depth=0):
        nonlocal findings_json
        if path is None:
            path = {"etherification": False, "halogenation": False, "reduction": False, "depth": {}}

        # Make a copy of the current path to avoid modifying parent paths
        current_path = path.copy()
        current_path["depth"] = path["depth"].copy()

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            print(f"Depth {depth} - Reaction: {rsmi}")

            # Check for etherification
            for r in ETHERIFICATION_REACTIONS:
                if checker.check_reaction(r, rsmi):
                    print(f"Found etherification at depth {depth}")
                    current_path["etherification"] = True
                    current_path["depth"]["etherification"] = depth
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)
                    break

            # Check for halogenation (OH to I/Br/Cl conversion)
            for r in ALCOHOL_HALOGENATION_REACTIONS:
                if checker.check_reaction(r, rsmi):
                    print(f"Found halogenation at depth {depth}")
                    current_path["halogenation"] = True
                    current_path["depth"]["halogenation"] = depth
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)
                    break

            # Check for carboxylic acid/ester reduction
            for r in ACID_ESTER_REDUCTION_REACTIONS:
                if checker.check_reaction(r, rsmi):
                    print(f"Found carboxylic acid/ester reduction at depth {depth}")
                    current_path["reduction"] = True
                    current_path["depth"]["reduction"] = depth
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)
                    break

        # If we have all three transformations, add this path to valid paths
        if (
            current_path["etherification"]
            and current_path["halogenation"]
            and current_path["reduction"]
        ):
            # Check if the transformations are in the correct order (by depth)
            eth_depth = current_path["depth"].get("etherification", float("inf"))
            hal_depth = current_path["depth"].get("halogenation", float("inf"))
            red_depth = current_path["depth"].get("reduction", float("inf"))

            # The correct order is: etherification (latest/lowest depth) -> halogenation -> reduction (earliest/highest depth)
            if eth_depth < hal_depth < red_depth:
                print(
                    f"Found complete path with correct order: etherification({eth_depth}) -> halogenation({hal_depth}) -> reduction({red_depth})"
                )
                valid_paths.append(current_path)

        # Process children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, current_path, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if any valid paths were found
    strategy_present = len(valid_paths) > 0

    if strategy_present:
        print("Found at least one valid path with the correct sequence of transformations")
        # Add structural constraints if the strategy is present
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "etherification_step",
                    "halogenation_step",
                    "reduction_step"
                ],
                "description": "A valid synthetic route must contain at least one reaction from each of the three categories: etherification, alcohol halogenation, and acid/ester reduction."
            }
        })
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_events": [
                    "reduction_step",
                    "halogenation_step",
                    "etherification_step"
                ],
                "description": "The reaction steps must occur in a specific order from earliest to latest: reduction, then halogenation, then etherification."
            }
        })
    else:
        print("No valid paths found with the correct sequence of transformations")

    return strategy_present, findings_json
