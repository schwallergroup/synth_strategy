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

def main(route) -> Tuple[bool, Dict]:
    """This function identifies syntheses where a morpholine ring is either
    formed and retained in the final product, or formed and subsequently
    opened in a later step. The check ensures the opening happens after the
    formation.
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

    # Track the stages of the synthesis
    morpholine_formation_stage = None
    morpholine_opening_stage = None

    # Track the presence of morpholine in the final product
    final_product_has_morpholine = False

    result = False # Initialize the main boolean result

    def dfs_traverse(node, depth=0):
        nonlocal morpholine_formation_stage, morpholine_opening_stage, final_product_has_morpholine, findings_json

        # Check if this is the final product (depth 0)
        if depth == 0 and node["type"] == "mol":
            if checker.check_ring("morpholine", node["smiles"]):
                final_product_has_morpholine = True
                findings_json["atomic_checks"]["ring_systems"].append("morpholine")

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for morpholine formation (ring closure)
                if checker.check_ring("morpholine", product) and not any(
                    checker.check_ring("morpholine", r) for r in reactants
                ):
                    if morpholine_formation_stage is None or depth < morpholine_formation_stage:
                        morpholine_formation_stage = depth
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                        if "morpholine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("morpholine")

                # Check for morpholine ring opening
                if any(
                    checker.check_ring("morpholine", r) for r in reactants
                ) and not checker.check_ring("morpholine", product):
                    if morpholine_opening_stage is None or depth < morpholine_opening_stage:
                        morpholine_opening_stage = depth
                        findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")
                        if "morpholine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("morpholine")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Check if key stages were found
    found_morpholine_formation = morpholine_formation_stage is not None
    found_morpholine_opening = morpholine_opening_stage is not None

    # Case 1: Final product has a morpholine that was formed during the synthesis.
    if final_product_has_morpholine:
        if found_morpholine_formation:
            result = True
            # Add structural constraint: positional (morpholine at last_stage)
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "morpholine", "position": "last_stage"}})
    # Case 2: Final product does not have a morpholine, but one was formed and later opened.
    # The opening must happen at a smaller depth (later in the synthesis) than the formation.
    else:
        # Add structural constraint: negation (morpholine not in final_product)
        if found_morpholine_formation or found_morpholine_opening: # Only add if there was some activity related to morpholine
            findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "morpholine", "scope": "final_product"}})

        if (
            found_morpholine_formation
            and found_morpholine_opening
            and morpholine_opening_stage < morpholine_formation_stage
        ):
            result = True
            # Add structural constraint: co-occurrence (ring_formation and ring_destruction)
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["ring_formation", "ring_destruction"], "scope": "morpholine"}})
            # Add structural constraint: sequence (ring_formation then ring_destruction)
            findings_json["structural_constraints"].append({"type": "sequence", "details": {"first": "ring_formation", "second": "ring_destruction", "target": "morpholine"}})

    return result, findings_json
