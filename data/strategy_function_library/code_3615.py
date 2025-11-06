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


AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann to ester",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis follows a convergent approach with a late-stage amide coupling.
    The function identifies an amide formation reaction occurring within the last two synthetic steps.
    It checks for specific named reactions from the AMIDE_COUPLING_REACTIONS list or, alternatively,
    for a reaction between a carboxylic acid/acyl halide and an amine. The strategy is confirmed
    if the synthesis tree has multiple branches and at least one of the reactants in the coupling
    step is complex (defined as >8 atoms).
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

    result = False

    # Track branch information for convergent pattern detection
    branch_info = {
        "has_multiple_branches": False,
        "branch_depths": [],
        "amide_coupling_depth": None,
        "complex_fragments_at_coupling": 0,
    }

    def dfs_traverse(node, depth=0, branch_id=0):
        nonlocal branch_info, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Only check for the coupling if one hasn't been found yet to prioritize the latest-stage one
            if branch_info["amide_coupling_depth"] is None and depth <= 2:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                is_amide_coupling = False
                for reaction_type in AMIDE_COUPLING_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_amide_coupling = True
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                if not is_amide_coupling and len(reactants) >= 2:
                    has_amide = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )
                    if has_amide:
                        if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"] and checker.check_fg("Primary amide", product):
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                        if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"] and checker.check_fg("Secondary amide", product):
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                        if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"] and checker.check_fg("Tertiary amide", product):
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                        amide_in_reactants = any(
                            checker.check_fg("Primary amide", r)
                            or checker.check_fg("Secondary amide", r)
                            or checker.check_fg("Tertiary amide", r)
                            for r in reactants
                        )
                        if amide_in_reactants:
                            # This corresponds to the 'negation' constraint
                            pass # We don't add it to findings_json if it's a negation
                        else:
                            if {"type": "negation", "details": {"target": "amide_group_in_reactants_of_amide_formation"}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "amide_group_in_reactants_of_amide_formation"}})

                        has_acid = False
                        has_amine = False
                        acid_reactant_idx = -1
                        amine_reactant_idx = -1
                        for i, r in enumerate(reactants):
                            if checker.check_fg("Carboxylic acid", r):
                                has_acid = True
                                acid_reactant_idx = i
                                if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                            if checker.check_fg("Acyl halide", r):
                                has_acid = True
                                acid_reactant_idx = i
                                if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")
                            if (
                                checker.check_fg("Primary amine", r)
                                or checker.check_fg("Secondary amine", r)
                                or checker.check_fg("Aniline", r)
                            ):
                                has_amine = True
                                amine_reactant_idx = i
                                if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"] and checker.check_fg("Primary amine", r):
                                    findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                                if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"] and checker.check_fg("Secondary amine", r):
                                    findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                                if "Aniline" not in findings_json["atomic_checks"]["functional_groups"] and checker.check_fg("Aniline", r):
                                    findings_json["atomic_checks"]["functional_groups"].append("Aniline")

                        if (
                            has_acid
                            and has_amine
                            and acid_reactant_idx != amine_reactant_idx
                            and not amide_in_reactants
                        ):
                            is_amide_coupling = True
                            if "amide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("amide_formation")

                if is_amide_coupling and len(reactants) >= 2:
                    branch_info["amide_coupling_depth"] = depth
                    if {"type": "positional", "details": {"target": "amide_coupling", "position": "depth <= 2"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amide_coupling", "position": "depth <= 2"}})

                    complex_fragments = 0
                    for reactant in reactants:
                        react_mol = Chem.MolFromSmiles(reactant)
                        if react_mol and react_mol.GetNumAtoms() > 8:
                            complex_fragments += 1
                    branch_info["complex_fragments_at_coupling"] = complex_fragments

        if len(node.get("children", [])) > 1:
            branch_info["has_multiple_branches"] = True
            branch_info["branch_depths"].append(depth)

        for i, child in enumerate(node.get("children", [])):
            new_branch_id = branch_id
            if len(node.get("children", [])) > 1:
                new_branch_id = branch_id + i + 1
            
            # New depth calculation logic
            if node["type"] == "reaction":
                next_depth = depth
            else: # Assuming 'chemical' or other types
                next_depth = depth + 1
            
            dfs_traverse(child, next_depth, new_branch_id)

    dfs_traverse(route)

    if (
        branch_info["has_multiple_branches"]
        and branch_info["complex_fragments_at_coupling"] >= 1
        and branch_info["amide_coupling_depth"] is not None
        and branch_info["amide_coupling_depth"] <= 2
    ):
        result = True
        if branch_info["has_multiple_branches"]:
            if {"type": "co-occurrence", "details": {"targets": ["convergent_synthesis", "late_stage_amide_coupling"]}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["convergent_synthesis", "late_stage_amide_coupling"]}})
        if branch_info["complex_fragments_at_coupling"] >= 1:
            if {"type": "count", "details": {"target": "complex_reactant_in_amide_coupling", "operator": ">=", "value": 1}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "count", "details": {"target": "complex_reactant_in_amide_coupling", "operator": ">=", "value": 1}})

    return result, findings_json
