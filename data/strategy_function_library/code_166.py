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


# Refactored lists of specific reaction names
BENZIMIDAZOLE_FORMATION_REACTIONS = [
    "benzimidazole_formation from aldehyde",
    "benzimidazole_formation from acyl halide",
    "benzimidazole_formation from ester/carboxylic acid",
    "{benzimidazole_derivatives_carboxylic-acid/ester}",
    "{benzimidazole_derivatives_aldehyde}",
]

BOC_PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
]

BOC_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a multi-part strategy defined by: (1) A late-stage (final two steps) formation of a benzimidazole ring, identified by a list of specific reaction names. (2) One or more Boc protection or deprotection reactions occurring in the early stages, also identified by specific reaction name lists. This function relies on the module-level constants BENZIMIDAZOLE_FORMATION_REACTIONS, BOC_PROTECTION_REACTIONS, and BOC_DEPROTECTION_REACTIONS.
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

    # Track key features
    has_late_benzimidazole_coupling = False
    has_boc_protection = False
    has_nitrile_intermediate = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_benzimidazole_coupling, has_boc_protection, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for benzimidazole formation at late stage (depth 0-1)
            if depth <= 1:
                for r in BENZIMIDAZOLE_FORMATION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        has_late_benzimidazole_coupling = True
                        if r not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)
                        # Add structural constraint for late-stage benzimidazole
                        if {"type": "positional", "details": {"target": BENZIMIDAZOLE_FORMATION_REACTIONS, "position": "last_two_stages"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": BENZIMIDAZOLE_FORMATION_REACTIONS, "position": "last_two_stages"}})
                        break

            # Check for Boc protection/deprotection in early stages (depth > 1)
            if depth > 1:
                boc_found = False
                for r in BOC_PROTECTION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        has_boc_protection = True
                        boc_found = True
                        if r not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)
                        break
                if not boc_found:
                    for r in BOC_DEPROTECTION_REACTIONS:
                        if checker.check_reaction(r, rsmi):
                            has_boc_protection = True
                            boc_found = True
                            if r not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(r)
                            break
                # Add structural constraint for early-stage Boc
                if boc_found and {"type": "positional", "details": {"target": BOC_PROTECTION_REACTIONS + BOC_DEPROTECTION_REACTIONS, "position": "early_stages"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": BOC_PROTECTION_REACTIONS + BOC_DEPROTECTION_REACTIONS, "position": "early_stages"}})

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have late-stage heterocycle coupling and earlier protection
    strategy_present = has_late_benzimidazole_coupling and has_boc_protection

    # Additional confidence if we also have nitrile intermediate
    if strategy_present and has_nitrile_intermediate:
        print(
            "Confirmed full late-stage heterocycle coupling strategy with protection and nitrile intermediate"
        )
    elif strategy_present:
        print("Detected late-stage heterocycle coupling strategy with protection")
    else:
        if has_late_benzimidazole_coupling:
            print("Found late-stage benzimidazole coupling but no Boc protection")
        if has_boc_protection:
            print("Found Boc protection but no late-stage benzimidazole coupling")

    return strategy_present, findings_json
