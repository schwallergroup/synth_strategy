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


BOC_PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
]

NITROGEN_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Tert-butyl deprotection of amine",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "N-glutarimide deprotection",
    "Phthalimide deprotection",
]

UREA_FORMATION_REACTIONS = [
    "Urea synthesis via isocyanate and primary amine",
    "Urea synthesis via isocyanate and secondary amine",
    "Urea synthesis via isocyanate and diazo",
    "Urea synthesis via isocyanate and sulfonamide",
    "urea",
    "{urea}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects two strategies for urea synthesis: 1) A protection-deprotection sequence where a nitrogen is protected (specifically with Boc), later deprotected, and finally converted to a urea. 2) Direct, late-stage formation of a urea. The function checks for specific, known reaction types for each step, as defined in the `BOC_PROTECTION_REACTIONS`, `NITROGEN_DEPROTECTION_REACTIONS`, and `UREA_FORMATION_REACTIONS` lists. The strategy is only flagged if a urea group is present in the final product.
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

    # Track key events in the synthesis
    protection_depth = -1
    deprotection_depth = -1
    urea_formation_depth = -1
    final_product_has_urea = False

    def dfs_traverse(node, depth=0):
        nonlocal protection_depth, deprotection_depth, urea_formation_depth, final_product_has_urea, findings_json

        # Check molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if the final product (depth 0) has a urea group
            if depth == 0:
                if checker.check_fg("Urea", mol_smiles):
                    final_product_has_urea = True
                    findings_json["atomic_checks"]["functional_groups"].append("Urea")
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Urea", "target_type": "functional_group", "position": "final_product"}})
                    print(f"Final product has urea: {mol_smiles}")

        # Check reaction nodes
        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")

            # Check for protection reactions (Boc or other common protecting groups)
            for rxn in BOC_PROTECTION_REACTIONS:
                if checker.check_reaction(rxn, rsmi):
                    if protection_depth == -1 or depth > protection_depth:
                        protection_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        print(f"Protection at depth {depth}: {rsmi}")
                    break # Only need to find one matching reaction

            # Check for deprotection reactions
            for rxn in NITROGEN_DEPROTECTION_REACTIONS:
                if checker.check_reaction(rxn, rsmi):
                    if deprotection_depth == -1 or depth < deprotection_depth:
                        deprotection_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        print(f"Deprotection at depth {depth}: {rsmi}")
                    break # Only need to find one matching reaction

            # Check for urea formation reactions
            if not any(checker.check_fg("Urea", r) for r in reactants) and checker.check_fg(
                "Urea", product_part
            ):
                # Check specific urea formation reactions
                for rxn in UREA_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        if urea_formation_depth == -1 or depth < urea_formation_depth:
                            urea_formation_depth = depth
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            print(f"Urea formation at depth {depth}: {rsmi}")
                        break # Only need to find one matching reaction

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is 'mol' (chemical), depth increases
                new_depth = depth + 1
            # If current node is 'reaction', depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Determine if the strategy is detected
    has_protection = protection_depth != -1
    has_deprotection = deprotection_depth != -1
    has_urea_formation = urea_formation_depth != -1

    # Check for protection-deprotection-urea sequence
    protection_sequence = (
        has_protection
        and has_deprotection
        and has_urea_formation
        and protection_depth > deprotection_depth
        and deprotection_depth > urea_formation_depth
    )
    if protection_sequence:
        findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["Boc amine protection", "Nitrogen deprotection", "urea_formation_de_novo"]}})

    # Check for direct urea formation (without protection/deprotection)
    direct_urea_formation = has_urea_formation and urea_formation_depth <= 2
    if direct_urea_formation:
        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "urea_formation_de_novo", "target_type": "reaction", "position_constraint": "depth <= 2"}})

    # Strategy is detected if either condition is met and final product has urea
    strategy_detected = (protection_sequence or direct_urea_formation) and final_product_has_urea

    print(f"Protection detected: {has_protection} at depth {protection_depth}")
    print(f"Deprotection detected: {has_deprotection} at depth {deprotection_depth}")
    print(f"Urea formation detected: {has_urea_formation} at depth {urea_formation_depth}")
    print(f"Final product has urea: {final_product_has_urea}")
    print(f"Protection sequence: {protection_sequence}")
    print(f"Direct urea formation: {direct_urea_formation}")
    print(f"Protection-replacement urea strategy: {strategy_detected}")

    # Remove duplicate entries from lists in findings_json
    for key in findings_json["atomic_checks"]:
        findings_json["atomic_checks"][key] = list(set(findings_json["atomic_checks"][key]))
    
    # For structural constraints, we need to ensure uniqueness based on content
    unique_constraints = []
    seen_constraints = set()
    for constraint in findings_json["structural_constraints"]:
        constraint_str = str(constraint) # Convert dict to string for hashing
        if constraint_str not in seen_constraints:
            unique_constraints.append(constraint)
            seen_constraints.add(constraint_str)
    findings_json["structural_constraints"] = unique_constraints

    return strategy_detected, findings_json
