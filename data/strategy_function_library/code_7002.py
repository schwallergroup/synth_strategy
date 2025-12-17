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


STRATEGIC_RING_FORMING_REACTIONS = [
    "Diels-Alder",
    "Huisgen",
    "Paal-Knorr pyrrole",
    "Pictet-Spengler",
    "Fischer indole",
    "Friedlaender chinoline",
    "benzofuran",
    "benzothiophene",
    "indole",
    "oxadiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy where fragments are combined early in the synthesis,
    followed by ring formation in the middle stage.

    In retrosynthetic analysis:
    - Early stage = high depth values
    - Late stage = low depth values
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

    # Initialize tracking variables
    fragment_combination_depth = -1
    ring_formation_depth = -1
    max_depth = -1

    # Track fragment combination reactions for later verification
    fragment_combinations = []

    def dfs_traverse(node, depth=0):
        nonlocal fragment_combination_depth, ring_formation_depth, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for fragment combination (more than one reactant)
            if len(reactants_smiles) > 1 and all(reactants_smiles):
                # In retrosynthetic analysis, this is a disconnection of a single molecule into fragments
                # Store the depth and the reaction for later verification
                if fragment_combination_depth == -1 or depth > fragment_combination_depth:
                    fragment_combination_depth = depth
                    fragment_combinations.append((depth, rsmi))
                    if "fragment_combination" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("fragment_combination")
                print(f"Fragment combination detected at depth {depth}")

            # Check for ring formation
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and all(reactants):
                # Get the number of rings in product and reactants
                product_ring_count = product.GetRingInfo().NumRings()
                reactant_ring_counts = [
                    r.GetRingInfo().NumRings() for r in reactants if r is not None
                ]
                max_reactant_rings = max(reactant_ring_counts, default=0)

                # In retrosynthetic analysis, this is breaking a ring in the product
                if product_ring_count > max_reactant_rings:
                    # Check if this is a specific, strategic ring-forming reaction
                    for rxn_name in STRATEGIC_RING_FORMING_REACTIONS:
                        if checker.check_reaction(rxn_name, rsmi):
                            if ring_formation_depth == -1 or depth < ring_formation_depth:
                                ring_formation_depth = depth
                            if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                            print(
                                f"Ring formation detected at depth {depth}, reaction type: {rxn_name}"
                            )
                            break

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    strategy_detected = False
    # Check if the strategy is present
    if max_depth > 0:
        early_third = max_depth // 3
        middle_third = max_depth * 2 // 3

        # In retrosynthetic analysis:
        # - Early stage = high depth values (max_depth - early_third to max_depth)
        # - Middle stage = middle depth values (max_depth - middle_third to max_depth - early_third)
        # - Late stage = low depth values (0 to max_depth - middle_third)

        # Fragment combination should be in early stage (high depth)
        early_fragment = (
            fragment_combination_depth >= (max_depth - early_third)
            if fragment_combination_depth != -1
            else False
        )
        if early_fragment:
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "fragment_combination",
                    "position": "early_stage"
                }
            })

        # Ring formation should be in middle stage (middle depth)
        middle_ring = (
            (
                ring_formation_depth < (max_depth - early_third)
                and ring_formation_depth >= (max_depth - middle_third)
            )
            if ring_formation_depth != -1
            else False
        )
        if middle_ring:
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "ring_formation",
                    "position": "middle_stage"
                }
            })

        # Verify that ring formation happens after fragment combination
        correct_sequence = (
            ring_formation_depth < fragment_combination_depth
            if (ring_formation_depth != -1 and fragment_combination_depth != -1)
            else False
        )
        if correct_sequence:
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "fragment_combination",
                    "after": "ring_formation"
                }
            })

        strategy_detected = early_fragment and middle_ring and correct_sequence

        if strategy_detected:
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "fragment_combination",
                        "ring_formation"
                    ]
                }
            })
            print("Early fragment combination with middle-stage ring formation strategy detected")
            print(
                f"Fragment combination depth: {fragment_combination_depth}, Ring formation depth: {ring_formation_depth}"
            )
            print(
                f"Early threshold: {max_depth - early_third}, Middle threshold: {max_depth - middle_third}, Max depth: {max_depth}"
            )
        else:
            print("Strategy not detected:")
            print(
                f"Early fragment: {early_fragment}, Middle ring: {middle_ring}, Correct sequence: {correct_sequence}"
            )
            print(
                f"Fragment combination depth: {fragment_combination_depth}, Ring formation depth: {ring_formation_depth}"
            )
            print(
                f"Early threshold: {max_depth - early_third}, Middle threshold: {max_depth - middle_third}, Max depth: {max_depth}"
            )

        return strategy_detected, findings_json

    return False, findings_json