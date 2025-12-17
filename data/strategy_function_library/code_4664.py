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


AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with ammonia to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Carboxylic acid to amide conversion",
    "Aminolysis of esters",
    "Nitrile to amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage amide formation as the final synthetic step. This is identified either by matching a predefined list of named amide formation reactions (e.g., Schotten-Baumann, aminolysis of esters, nitrile hydrolysis) or by a fallback analysis confirming the presence of amine (primary/secondary) and carboxylic acid derivative reactants leading to an amide product.
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

    final_step_is_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_amide_formation, findings_json

        # Check if this is the final reaction step (depth=1) and it's a reaction
        if depth == 1 and node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                try:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    reactants_smiles = rsmi.split(">")[0]
                    product_smiles = rsmi.split(">")[-1]

                    # Check if product contains an amide group
                    has_product_amide = False
                    if checker.check_fg("Primary amide", product_smiles):
                        has_product_amide = True
                        if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                    if checker.check_fg("Secondary amide", product_smiles):
                        has_product_amide = True
                        if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                    if checker.check_fg("Tertiary amide", product_smiles):
                        has_product_amide = True
                        if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                    if not has_product_amide:
                        return

                    # Check if any of the amide formation reactions match
                    for reaction_type in AMIDE_FORMATION_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            final_step_is_amide_formation = True
                            if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            # Add positional constraint if a named reaction matches
                            if {"type": "positional", "details": {"target": "amide_formation", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amide_formation", "position": "last_stage"}})
                            return

                    # If no specific reaction matched, check for functional group changes
                    reactants_list = reactants_smiles.split(".")

                    # Check for carboxylic acid derivatives in reactants
                    has_acid = False
                    if any(checker.check_fg("Carboxylic acid", r) for r in reactants_list):
                        has_acid = True
                        if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                    has_acyl_halide = False
                    if any(checker.check_fg("Acyl halide", r) for r in reactants_list):
                        has_acyl_halide = True
                        if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")

                    has_ester = False
                    if any(checker.check_fg("Ester", r) for r in reactants_list):
                        has_ester = True
                        if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ester")

                    has_anhydride = False
                    if any(checker.check_fg("Anhydride", r) for r in reactants_list):
                        has_anhydride = True
                        if "Anhydride" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Anhydride")

                    # Check for amines in reactants that can form amides
                    has_amine = False
                    if any(checker.check_fg("Primary amine", r) for r in reactants_list):
                        has_amine = True
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    if any(checker.check_fg("Secondary amine", r) for r in reactants_list):
                        has_amine = True
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                    if any(checker.check_fg("Aniline", r) for r in reactants_list):
                        has_amine = True
                        if "Aniline" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aniline")

                    # Check if we have the necessary components for amide formation
                    if (
                        has_amine
                        and (
                            has_acid or has_acyl_halide or has_ester or has_anhydride
                        )
                    ) and has_product_amide:
                        final_step_is_amide_formation = True
                        # Add co-occurrence constraint if fallback logic matches
                        if {"type": "co-occurrence", "details": {"description": "As a fallback, checks for the co-occurrence of specific functional group categories in the reactants and products of a single reaction step. An amine-type reactant and an acid-derivative-type reactant must lead to an amide-type product.", "targets": ["amine_reactant", "acid_derivative_reactant", "amide_product"]}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"description": "As a fallback, checks for the co-occurrence of specific functional group categories in the reactants and products of a single reaction step. An amine-type reactant and an acid-derivative-type reactant must lead to an amide-type product.", "targets": ["amine_reactant", "acid_derivative_reactant", "amide_product"]}})
                        # Add positional constraint if fallback logic matches
                        if {"type": "positional", "details": {"target": "amide_formation", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amide_formation", "position": "last_stage"}})

                except Exception as e:
                    print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return final_step_is_amide_formation, findings_json
