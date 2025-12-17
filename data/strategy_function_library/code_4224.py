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


ESTERIFICATION_REACTION_NAMES = [
    "Esterification of Carboxylic Acids",
    "Schotten-Baumann to ester",
    "O-alkylation of carboxylic acids with diazo compounds",
    "Transesterification",
    "Acetic anhydride and alcohol to ester",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage (depth <= 1) ester formation. It checks for a specific list of named esterification reactions: 'Esterification of Carboxylic Acids', 'Schotten-Baumann to ester', 'O-alkylation of carboxylic acids with diazo compounds', 'Transesterification', and 'Acetic anhydride and alcohol to ester'. It also identifies general esterifications by checking for an alcohol and a suitable acyl donor (carboxylic acid, anhydride, acyl halide) in the reactants, with confirmation of an ester in the products.
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

    has_late_stage_esterification = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_esterification, findings_json

        if node["type"] == "reaction" and depth <= 1:
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Direct check for named esterification reactions
                is_named_esterification = False
                for name in ESTERIFICATION_REACTION_NAMES:
                    if checker.check_reaction(name, rsmi):
                        is_named_esterification = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)

                # Check for reactant functional groups for a general esterification
                has_alcohol = False
                for r in reactants:
                    if not r: continue
                    if checker.check_fg("Primary alcohol", r):
                        has_alcohol = True
                        if "Primary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                    if checker.check_fg("Secondary alcohol", r):
                        has_alcohol = True
                        if "Secondary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
                    if checker.check_fg("Tertiary alcohol", r):
                        has_alcohol = True
                        if "Tertiary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")
                    if checker.check_fg("Aromatic alcohol", r):
                        has_alcohol = True
                        if "Aromatic alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aromatic alcohol")

                has_carboxyl = False
                for r in reactants:
                    if not r: continue
                    if checker.check_fg("Carboxylic acid", r):
                        has_carboxyl = True
                        if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                has_anhydride = False
                for r in reactants:
                    if not r: continue
                    if checker.check_fg("Anhydride", r):
                        has_anhydride = True
                        if "Anhydride" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Anhydride")

                has_acyl_halide = False
                for r in reactants:
                    if not r: continue
                    if checker.check_fg("Acyl halide", r):
                        has_acyl_halide = True
                        if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")

                # Confirm ester product for the general case
                has_ester_product = False
                if product and checker.check_fg("Ester", product):
                    has_ester_product = True
                    if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ester")

                is_general_esterification = (
                    (
                        (has_alcohol and has_carboxyl)
                        or (has_alcohol and has_anhydride)
                        or (has_alcohol and has_acyl_halide)
                    )
                    and has_ester_product
                )

                if is_named_esterification or is_general_esterification:
                    has_late_stage_esterification = True
                    # Add the structural constraint if the condition is met
                    if {"type": "positional", "details": {"target": "ester_formation", "position": "last_two_stages"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ester_formation", "position": "last_two_stages"}})

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return has_late_stage_esterification, findings_json
