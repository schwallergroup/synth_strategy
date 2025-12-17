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


ARYL_ETHER_REACTIONS = [
    "Williamson Ether Synthesis",
    "Mitsunobu aryl ether",
    "Chan-Lam etherification",
    "Ullmann-Goldberg Substitution aryl alcohol",
]

AZIDE_PRECURSOR_FGS = [
    "Primary halide",
    "Secondary halide",
    "Tertiary halide",
    "Aromatic halide",
    "Triflate",
    "Mesylate",
    "Tosylate",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route involves all three key features:
    1.  Late-stage azide introduction
    2.  Epoxide ring opening
    3.  Early aryl ether formation
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

    azide_introduced = False
    azide_depth = -1
    epoxide_opening_detected = False
    epoxide_depth = -1
    aryl_ether_formed = False
    aryl_ether_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal azide_introduced, azide_depth, epoxide_opening_detected, epoxide_depth, aryl_ether_formed, aryl_ether_depth, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check for azide introduction
            product_has_azide = checker.check_fg("Azide", product)
            if product_has_azide:
                findings_json["atomic_checks"]["functional_groups"].append("Azide")

            reactants_have_azide = any(
                checker.check_fg("Azide", reactant) for reactant in reactants if reactant
            )

            if product_has_azide and not reactants_have_azide:
                print(f"Azide introduced at depth {depth}, reaction: {rsmi}")
                azide_introduced = True
                azide_depth = depth

            # Check for epoxide ring opening
            reactants_have_epoxide = False
            for reactant in reactants:
                if reactant and checker.check_ring("oxirane", reactant):
                    reactants_have_epoxide = True
                    findings_json["atomic_checks"]["ring_systems"].append("oxirane")
                    break

            product_has_epoxide = checker.check_ring("oxirane", product)

            if reactants_have_epoxide and not product_has_epoxide:
                print(f"Epoxide ring opening detected at depth {depth}, reaction: {rsmi}")
                epoxide_opening_detected = True
                epoxide_depth = depth

            # Check for aryl ether formation via specific named reactions
            for rxn in ARYL_ETHER_REACTIONS:
                if checker.check_reaction(rxn, rsmi):
                    print(f"Aryl ether formation via specific reaction at depth {depth}")
                    aryl_ether_formed = True
                    aryl_ether_depth = depth
                    findings_json["atomic_checks"]["named_reactions"].append(rxn)
                    break

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(
        f"Summary - Azide: {azide_introduced} (depth {azide_depth}), Epoxide: {epoxide_opening_detected} (depth {epoxide_depth}), Aryl ether: {aryl_ether_formed} (depth {aryl_ether_depth})"
    )

    # Check if all three features are present and in the correct order
    # Late-stage is defined as depth <= 2, early-stage as depth >= 3
    late_stage_azide = azide_introduced and azide_depth <= 2
    mid_stage_epoxide = epoxide_opening_detected and (epoxide_depth >= 1 and epoxide_depth <= 4)
    early_aryl_ether = aryl_ether_formed and aryl_ether_depth >= 3

    print(
        f"Stage checks - Late azide: {late_stage_azide}, Mid epoxide: {mid_stage_epoxide}, Early aryl ether: {early_aryl_ether}"
    )

    result = False

    # Check for the combined strategy
    if late_stage_azide and mid_stage_epoxide and early_aryl_ether:
        print(
            "Combined strategy detected: Late-stage azide introduction, epoxide ring-opening, with early aryl ether formation"
        )
        result = True
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"description": "The primary strategy requires the co-occurrence of three key events: aryl ether formation, epoxide ring opening, and azide introduction.", "targets": ["aryl_ether_formation", "epoxide_ring_opening", "azide_introduction"]}})
        findings_json["structural_constraints"].append({"type": "positional", "details": {"description": "Aryl ether formation must occur in an early stage of the synthesis.", "target": "aryl_ether_formation", "position": "early_stage (depth >= 3)"}})
        findings_json["structural_constraints"].append({"type": "positional", "details": {"description": "Epoxide ring opening should occur in the middle stages of the synthesis.", "target": "epoxide_ring_opening", "position": "mid_stage (depth 1-4)"}})
        findings_json["structural_constraints"].append({"type": "positional", "details": {"description": "Azide introduction must be a late-stage transformation.", "target": "azide_introduction", "position": "late_stage (depth <= 2)"}})

    # If we're missing the azide introduction but have the other features, check for potential azide precursors
    if not azide_introduced and mid_stage_epoxide and early_aryl_ether:
        # Check if there's a halide or other good leaving group in the final product that could be converted to azide
        target_mol = route["smiles"]
        has_leaving_group = False
        for fg in AZIDE_PRECURSOR_FGS:
            if checker.check_fg(fg, target_mol):
                has_leaving_group = True
                findings_json["atomic_checks"]["functional_groups"].append(fg)
                break

        print(f"Target molecule has leaving group: {has_leaving_group}")

        if has_leaving_group:
            print(
                "Combined strategy detected with potential for azide introduction in the final step"
            )
            result = True
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"description": "The primary strategy requires the co-occurrence of three key events: aryl ether formation, epoxide ring opening, and azide introduction.", "targets": ["aryl_ether_formation", "epoxide_ring_opening", "azide_introduction"]}})
            findings_json["structural_constraints"].append({"type": "positional", "details": {"description": "Aryl ether formation must occur in an early stage of the synthesis.", "target": "aryl_ether_formation", "position": "early_stage (depth >= 3)"}})
            findings_json["structural_constraints"].append({"type": "positional", "details": {"description": "Epoxide ring opening should occur in the middle stages of the synthesis.", "target": "epoxide_ring_opening", "position": "mid_stage (depth 1-4)"}})
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"description": "An alternative strategy is accepted if the azide is not introduced, but an azide precursor functional group is present on the final product, along with the other key events.", "targets": ["aryl_ether_formation", "epoxide_ring_opening", "azide_precursor_on_final_product"]}})
            findings_json["structural_constraints"].append({"type": "negation", "details": {"description": "The alternative strategy is only valid if the azide introduction event has not occurred in the route.", "target": "azide_introduction"}})

    return result, findings_json