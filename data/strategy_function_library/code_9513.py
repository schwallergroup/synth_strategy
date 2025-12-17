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
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Schotten-Baumann_amide",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
]

ESTER_HYDROLYSIS_REACTIONS = [
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "COOH ethyl deprotection",
    "Deprotection of carboxylic acid",
]

AMINE_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Tert-butyl deprotection of amine",
    "Phthalimide deprotection",
    "N-glutarimide deprotection",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis featuring a late-stage amide coupling where at least two fragments
    undergo a deprotection/hydrolysis step prior to the final coupling. The specific reactions checked
    for amide coupling are defined in `AMIDE_COUPLING_REACTIONS`. The deprotection reactions checked
    are defined in `ESTER_HYDROLYSIS_REACTIONS` and `AMINE_DEPROTECTION_REACTIONS`.
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

    amide_coupling_node = None
    branches = {}
    result = False

    def dfs_traverse(node, depth=0, branch_id=None):
        nonlocal amide_coupling_node, result, findings_json

        if node.get("type") == "reaction":
            if "metadata" in node and "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]
                reactants = reactants_str.split(".")

                # Check for amide coupling at depth <= 1 (late stage)
                if depth <= 1 and not amide_coupling_node:
                    is_amide_coupling = False

                    for reaction_name in AMIDE_COUPLING_REACTIONS:
                        if checker.check_reaction(reaction_name, rsmi):
                            is_amide_coupling = True
                            if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                            break

                    if not is_amide_coupling:
                        has_acid = False
                        has_amine = False
                        for reactant in reactants:
                            if checker.check_fg("Carboxylic acid", reactant):
                                has_acid = True
                                if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                            elif checker.check_fg("Acyl halide", reactant):
                                has_acid = True
                                if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")
                            elif checker.check_fg("Primary amine", reactant):
                                has_amine = True
                                if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                            elif checker.check_fg("Secondary amine", reactant):
                                has_amine = True
                                if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

                        has_amide = False
                        if (
                            checker.check_fg("Primary amide", product_str)
                        ):
                            has_amide = True
                            if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                        if (
                            checker.check_fg("Secondary amide", product_str)
                        ):
                            has_amide = True
                            if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                        if (
                            checker.check_fg("Tertiary amide", product_str)
                        ):
                            has_amide = True
                            if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                        if has_acid and has_amine and has_amide:
                            is_amide_coupling = True

                    if is_amide_coupling and len(reactants) >= 2:
                        amide_coupling_node = node
                        if {"type": "positional", "details": {"target": "Amide coupling", "position": "late_stage", "min_reactants": 2}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Amide coupling", "position": "late_stage", "min_reactants": 2}})
                        for i, reactant in enumerate(reactants):
                            branch_id = f"branch_{i}"
                            branches[branch_id] = {
                                "has_deprotection": False,
                                "reactant": reactant,
                                "depth": depth,
                            }

                # Check for deprotection/hydrolysis in any branch
                if branch_id and branch_id in branches:
                    # Check for ester hydrolysis
                    for reaction_name in ESTER_HYDROLYSIS_REACTIONS:
                        if checker.check_reaction(reaction_name, rsmi):
                            branches[branch_id]["has_deprotection"] = True
                            if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                            break

                    # Check for amine deprotection
                    if not branches[branch_id]["has_deprotection"]:
                        for reaction_name in AMINE_DEPROTECTION_REACTIONS:
                            if checker.check_reaction(reaction_name, rsmi):
                                branches[branch_id]["has_deprotection"] = True
                                if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                                break

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node.get("type") != "reaction": # This means it's a 'chemical' node
            next_depth = depth + 1

        # Process children
        if node == amide_coupling_node and node.get("children"):
            for i, child in enumerate(node.get("children", [])):
                if i < len(branches):
                    child_branch_id = f"branch_{i}"
                    dfs_traverse(child, next_depth, child_branch_id)
                else:
                    dfs_traverse(child, next_depth, None)
        else:
            for child in node.get("children", []):
                dfs_traverse(child, next_depth, branch_id)

    dfs_traverse(route)

    if amide_coupling_node:
        deprotection_count = sum(1 for branch in branches.values() if branch["has_deprotection"])

        # A true convergent deprotection-coupling requires at least two branches to have deprotection steps.
        if deprotection_count >= 2:
            result = True
            if {"type": "count", "details": {"target": "precursor_branch_with_deprotection", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "count", "details": {"target": "precursor_branch_with_deprotection", "operator": ">=", "value": 2}})

    return result, findings_json
