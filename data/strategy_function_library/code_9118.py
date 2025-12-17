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


AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Schotten-Baumann to ester",
    "{Schotten-Baumann_amide}",
    "Acyl chloride with secondary amine to amide",
    "Ester with secondary amine to amide",
    "Acyl chloride with ammonia to amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthesis strategy that combines a late-stage amide coupling with a nitro reduction.
    It verifies that:
    1. An amide bond is formed in the final synthetic step (depth=1), identified by checking against a list of named reactions (AMIDE_COUPLING_REACTIONS) or by confirming the transformation of a carboxylic acid derivative and an amine into an amide.
    2. A nitro group is reduced to an amine at any point in the synthesis.
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

    has_amide_coupling_final_step = False
    has_nitro_reduction = False

    # This part is now obsolete due to the logic fix in dfs_traverse, but cannot be removed per instructions.
    final_step_depth = None
    def find_final_step(node, depth=0):
        nonlocal final_step_depth
        if node["type"] == "reaction":
            if final_step_depth is None or depth < final_step_depth:
                final_step_depth = depth
        for child in node.get("children", []):
            find_final_step(child, depth + 1)
    find_final_step(route)

    def dfs_traverse(node, depth=0, path=None):
        nonlocal has_amide_coupling_final_step, has_nitro_reduction, findings_json

        if path is None:
            path = []

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check for amide coupling at the final step (depth == 1)
                if depth == 1:
                    is_amide_reaction = False
                    for rxn in AMIDE_COUPLING_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            is_amide_reaction = True
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            break

                    has_acid_derivative = False
                    has_amine = False
                    for reactant in reactants:
                        if checker.check_fg("Carboxylic acid", reactant):
                            has_acid_derivative = True
                            if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                        if checker.check_fg("Acyl halide", reactant):
                            has_acid_derivative = True
                            if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")
                        if checker.check_fg("Ester", reactant):
                            has_acid_derivative = True
                            if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Ester")
                        if checker.check_fg("Anhydride", reactant):
                            has_acid_derivative = True
                            if "Anhydride" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Anhydride")

                        if checker.check_fg("Primary amine", reactant):
                            has_amine = True
                            if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        if checker.check_fg("Secondary amine", reactant):
                            has_amine = True
                            if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                        if checker.check_fg("Aniline", reactant):
                            has_amine = True
                            if "Aniline" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Aniline")

                    has_amide_product = False
                    if checker.check_fg("Primary amide", product_part):
                        has_amide_product = True
                        if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                    if checker.check_fg("Secondary amide", product_part):
                        has_amide_product = True
                        if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                    if checker.check_fg("Tertiary amide", product_part):
                        has_amide_product = True
                        if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                    if (
                        not is_amide_reaction
                        and has_acid_derivative
                        and has_amine
                        and has_amide_product
                    ):
                        is_amide_reaction = True
                        # If it's an implicit amide formation, add a generic tag
                        if "amide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("amide_formation")

                    if is_amide_reaction and has_amide_product:
                        has_amide_coupling_final_step = True
                        if {"type": "positional", "details": {"target": "amide_coupling", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amide_coupling", "position": "last_stage"}})

                # Check for nitro reduction at any depth
                has_nitro_reactant = False
                for reactant in reactants:
                    if checker.check_fg("Nitro group", reactant):
                        has_nitro_reactant = True
                        if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

                if has_nitro_reactant:
                    has_amine_product = False
                    if checker.check_fg("Primary amine", product_part):
                        has_amine_product = True
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    if checker.check_fg("Aniline", product_part):
                        has_amine_product = True
                        if "Aniline" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aniline")

                    has_nitro_product = checker.check_fg("Nitro group", product_part)

                    if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                        has_nitro_reduction = True
                        if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                    elif has_amine_product and not has_nitro_product and has_nitro_reactant:
                        has_nitro_reduction = True
                        # If it's an implicit nitro reduction, add a generic tag
                        if "nitro_reduction" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("nitro_reduction")

            except Exception as e:
                pass

        current_path = path + [node]

        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                dfs_traverse(child, depth, current_path)
            else: # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, depth + 1, current_path)

    dfs_traverse(route)

    result = has_amide_coupling_final_step and has_nitro_reduction

    if result:
        # Add the co-occurrence constraint if both conditions are met
        if {"type": "co-occurrence", "details": {"targets": ["amide_coupling", "nitro_reduction"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["amide_coupling", "nitro_reduction"]}})

    return result, findings_json
