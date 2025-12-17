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


CLICK_REACTION_NAMES = [
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen_Cu-catalyzed_1,4-subst",
    "Huisgen_Ru-catalyzed_1,5_subst",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
    "Azide-nitrile click cycloaddition to tetrazole",
    "Azide-nitrile click cycloaddition to triazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects an intermolecular click chemistry step. This is defined by a reaction that matches a defined list of click reaction types (e.g., Huisgen cycloadditions) and where the core reacting functional groups (e.g., azide and alkyne/nitrile) are located on separate reactant molecules.
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

    convergent_click_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_click_detected, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0]

            # Check if this is a click reaction from the defined list
            is_click_reaction = False
            for name in CLICK_REACTION_NAMES:
                if checker.check_reaction(name, rsmi):
                    is_click_reaction = True
                    if name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(name)
            
            if is_click_reaction:
                reactants = reactants_smiles.split(".")

                if len(reactants) >= 2:
                    # Record structural constraint for reactant count
                    if {"type": "count", "details": {"target": "reactants_in_click_reaction", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reactants_in_click_reaction", "operator": ">=", "value": 2}})

                    azide_fragments = []
                    alkyne_fragments = []
                    nitrile_fragments = []

                    for i, r in enumerate(reactants):
                        if checker.check_fg("Azide", r):
                            azide_fragments.append(i)
                            if "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Azide")
                        if checker.check_fg("Alkyne", r):
                            alkyne_fragments.append(i)
                            if "Alkyne" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Alkyne")
                        if checker.check_fg("Nitrile", r):
                            nitrile_fragments.append(i)
                            if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

                    azide_alkyne_click = (
                        azide_fragments
                        and alkyne_fragments
                        and set(azide_fragments).isdisjoint(set(alkyne_fragments))
                    )

                    azide_nitrile_click = (
                        azide_fragments
                        and nitrile_fragments
                        and set(azide_fragments).isdisjoint(set(nitrile_fragments))
                    )

                    if azide_alkyne_click or azide_nitrile_click:
                        convergent_click_detected = True
                        # Record structural constraint for co-occurrence
                        co_occurrence_constraint = {"type": "co-occurrence", "details": {"targets": ["click_reaction", "Azide", "Alkyne", "Nitrile"]}}
                        if co_occurrence_constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(co_occurrence_constraint)

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return convergent_click_detected, findings_json
