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


LATE_STAGE_TARGET_FGS = ["Urea", "Thiourea"]
LATE_STAGE_FORMATION_REACTIONS = [
    "Urea synthesis via isocyanate and primary amine",
    "Urea synthesis via isocyanate and secondary amine",
    "Urea synthesis via isocyanate and diazo",
    "Urea synthesis via isocyanate and sulfonamide",
    "{urea}",
    "thiourea",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the late-stage formation of specific functional groups, including Urea and Thiourea.
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

    late_urea_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_urea_formation, findings_json

        if node["type"] == "reaction" and depth == 1:
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_has_target_fg = False
                for fg in LATE_STAGE_TARGET_FGS:
                    if checker.check_fg(fg, product):
                        product_has_target_fg = True
                        if fg not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg)
                        break

                if product_has_target_fg:
                    is_target_reaction = False
                    for rxn in LATE_STAGE_FORMATION_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            is_target_reaction = True
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            break

                    if is_target_reaction:
                        target_fg_in_reactants = False
                        for reactant in reactants:
                            for fg in LATE_STAGE_TARGET_FGS:
                                if checker.check_fg(fg, reactant):
                                    target_fg_in_reactants = True
                                    # Do not add to findings_json here, as it's a negative check (absence)
                                    break
                            if target_fg_in_reactants:
                                break

                        if not target_fg_in_reactants:
                            late_urea_formation = True
                            # Add the 'functional_group_formation' reaction if it's not already there
                            if "functional_group_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("functional_group_formation")
                            # Add the structural constraint for late-stage formation
                            constraint = {
                                "type": "positional",
                                "details": {
                                    "target": "functional_group_formation",
                                    "position": "last_stage"
                                }
                            }
                            if constraint not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(constraint)

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return late_urea_formation, findings_json
