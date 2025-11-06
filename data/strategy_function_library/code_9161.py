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


AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with ammonia to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
    "Carboxylic acid to amide conversion",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the route involves late-stage amide formation by checking for specific, named reactions.
    Late-stage is defined as occurring in the first half of the synthesis steps (where depth=1 is the final step).
    The function confirms amide formation by checking for the presence of a primary, secondary, or tertiary amide in the product that was not present in all reactants.
    The specific reaction types checked are defined in the `AMIDE_FORMATION_REACTIONS` list.
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

    amide_formation_depth = None
    max_depth = 0
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_depth, max_depth, result, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            is_amide_formation = False
            for reaction_type in AMIDE_FORMATION_REACTIONS:
                if checker.check_reaction(reaction_type, rsmi):
                    is_amide_formation = True
                    if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    break

            product_has_amide_group = (
                checker.check_fg("Primary amide", product)
                or checker.check_fg("Secondary amide", product)
                or checker.check_fg("Tertiary amide", product)
            )

            if checker.check_fg("Primary amide", product) and "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
            if checker.check_fg("Secondary amide", product) and "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
            if checker.check_fg("Tertiary amide", product) and "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

            # Verify that an amide group exists in the product
            if is_amide_formation and product_has_amide_group:
                # Add co-occurrence constraint
                if {"type": "co-occurrence", "details": {"targets": ["any_amide_formation_reaction", "product_has_amide_group"]}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["any_amide_formation_reaction", "product_has_amide_group"]}})

                # Check that amide wasn't already in all reactants
                all_reactants_contain_amide_group = all(
                    checker.check_fg("Primary amide", r)
                    or checker.check_fg("Secondary amide", r)
                    or checker.check_fg("Tertiary amide", r)
                    for r in reactants
                )
                if not all_reactants_contain_amide_group:
                    # Add negation constraint
                    if {"type": "negation", "details": {"target": "all_reactants_contain_amide_group"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "all_reactants_contain_amide_group"}})
                    amide_formation_depth = depth

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if amide_formation_depth is not None:
        is_late_stage = amide_formation_depth <= (max_depth / 2)
        if is_late_stage:
            result = True
            # Add positional constraint
            if {"type": "positional", "details": {"target": "verified_amide_formation", "position": "first_half"}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "verified_amide_formation", "position": "first_half"}})

    return result, findings_json
