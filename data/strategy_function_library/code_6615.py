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
    This function detects a synthetic strategy involving acetal protection and deprotection.
    It looks for acetal formation in early stages and deprotection in later stages.
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

    protection_reactions = []
    deprotection_reactions = []
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, findings_json
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check for acetal formation (protection)
            is_acetal_formation_aldehyde_ketone = checker.check_reaction(
                "Aldehyde or ketone acetalization", rsmi
            )
            is_acetal_formation_diol = checker.check_reaction("Diol acetalization", rsmi)

            if is_acetal_formation_aldehyde_ketone:
                if "Aldehyde or ketone acetalization" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Aldehyde or ketone acetalization")
            if is_acetal_formation_diol:
                if "Diol acetalization" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Diol acetalization")

            is_acetal_formation = is_acetal_formation_aldehyde_ketone or is_acetal_formation_diol

            # Verify ketone/aldehyde in reactants and acetal in product
            has_ketone_aldehyde_reactants = False
            for r in reactants:
                if checker.check_fg("Ketone", r):
                    if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                    has_ketone_aldehyde_reactants = True
                if checker.check_fg("Aldehyde", r):
                    if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                    has_ketone_aldehyde_reactants = True

            has_acetal_product = checker.check_fg("Acetal/Ketal", product)
            if has_acetal_product:
                if "Acetal/Ketal" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Acetal/Ketal")

            if is_acetal_formation and has_ketone_aldehyde_reactants and has_acetal_product:
                protection_reactions.append((depth, rsmi))
                print(f"Acetal protection found at depth {depth}")

            # Check for acetal deprotection
            is_acetal_deprotection_aldehyde = checker.check_reaction("Acetal hydrolysis to aldehyde", rsmi)
            is_ketal_deprotection_ketone = checker.check_reaction("Ketal hydrolysis to ketone", rsmi)
            is_acetal_deprotection_diol = checker.check_reaction("Acetal hydrolysis to diol", rsmi)

            if is_acetal_deprotection_aldehyde:
                if "Acetal hydrolysis to aldehyde" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Acetal hydrolysis to aldehyde")
            if is_ketal_deprotection_ketone:
                if "Ketal hydrolysis to ketone" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Ketal hydrolysis to ketone")
            if is_acetal_deprotection_diol:
                if "Acetal hydrolysis to diol" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Acetal hydrolysis to diol")

            is_acetal_deprotection = (
                is_acetal_deprotection_aldehyde
                or is_ketal_deprotection_ketone
                or is_acetal_deprotection_diol
            )

            # Verify acetal in reactants and ketone/aldehyde in product
            has_acetal_reactants = False
            for r in reactants:
                if checker.check_fg("Acetal/Ketal", r):
                    if "Acetal/Ketal" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Acetal/Ketal")
                    has_acetal_reactants = True

            has_ketone_aldehyde_product = False
            if checker.check_fg("Ketone", product):
                if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                has_ketone_aldehyde_product = True
            if checker.check_fg("Aldehyde", product):
                if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                has_ketone_aldehyde_product = True

            if is_acetal_deprotection and has_acetal_reactants and has_ketone_aldehyde_product:
                deprotection_reactions.append((depth, rsmi))
                print(f"Acetal deprotection found at depth {depth}")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Determine early and late stages based on the maximum depth
    early_stage_threshold = max_depth // 2

    # In retrosynthesis, higher depths are early stages and lower depths are late stages
    protection_in_early_stage = any(
        depth > early_stage_threshold for depth, _ in protection_reactions
    )
    deprotection_in_late_stage = any(
        depth <= early_stage_threshold for depth, _ in deprotection_reactions
    )

    # Alternative check: ensure protection happens before deprotection
    sequential_strategy = False
    if protection_reactions and deprotection_reactions:
        min_protection_depth = min(depth for depth, _ in protection_reactions)
        max_deprotection_depth = max(depth for depth, _ in deprotection_reactions)
        sequential_strategy = min_protection_depth > max_deprotection_depth

    # Strategy is found if either condition is met
    strategy_found = (
        protection_in_early_stage and deprotection_in_late_stage
    ) or sequential_strategy

    print(f"Max depth: {max_depth}")
    print(f"Early stage threshold: {early_stage_threshold}")
    print(f"Protection reactions: {protection_reactions}")
    print(f"Deprotection reactions: {deprotection_reactions}")
    print(f"Protection in early stage: {protection_in_early_stage}")
    print(f"Deprotection in late stage: {deprotection_in_late_stage}")
    print(f"Sequential strategy: {sequential_strategy}")
    print(f"Strategy found: {strategy_found}")

    # Populate structural constraints
    if protection_reactions and deprotection_reactions:
        if {"type": "co-occurrence", "details": {"targets": ["acetal_protection", "acetal_deprotection"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["acetal_protection", "acetal_deprotection"]}})

    if sequential_strategy:
        if {"type": "sequence", "details": {"before": "acetal_protection", "after": "acetal_deprotection"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "sequence", "details": {"before": "acetal_protection", "after": "acetal_deprotection"}})

    if protection_in_early_stage:
        if {"type": "positional", "details": {"target": "acetal_protection", "position": "early_stage"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "acetal_protection", "position": "early_stage"}})

    if deprotection_in_late_stage:
        if {"type": "positional", "details": {"target": "acetal_deprotection", "position": "late_stage"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "acetal_deprotection", "position": "late_stage"}})

    return strategy_found, findings_json
