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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route involves late-stage tetrazole formation from nitrile
    via [3+2] cycloaddition with azide.
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

    tetrazole_formed = False
    nitrile_precursor = False
    azide_precursor = False
    tetrazole_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal tetrazole_formed, nitrile_precursor, azide_precursor, tetrazole_depth, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this reaction forms a tetrazole via azide-nitrile cycloaddition
            if checker.check_reaction("Azide-nitrile click cycloaddition to tetrazole", rsmi):
                tetrazole_formed = True
                tetrazole_depth = min(tetrazole_depth, depth)
                # The named reaction implies the presence of both precursors for this step.
                nitrile_precursor = True
                azide_precursor = True
                findings_json["atomic_checks"]["named_reactions"].append("Azide-nitrile click cycloaddition to tetrazole")
                # Since this reaction implies both precursors, we can add them here
                if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                if "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Azide")
                if "tetrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("tetrazole")

            # Alternative check: product has tetrazole and reactants have nitrile/azide
            elif checker.check_ring("tetrazole", product):
                if "tetrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("tetrazole")

                has_nitrile = False
                has_azide = False
                for reactant in reactants:
                    if checker.check_fg("Nitrile", reactant):
                        has_nitrile = True
                        if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                    if checker.check_fg("Azide", reactant):
                        has_azide = True
                        if "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Azide")

                # If we have both nitrile and azide, it's a strong indicator of the desired reaction.
                if has_nitrile and has_azide:
                    tetrazole_formed = True
                    tetrazole_depth = min(tetrazole_depth, depth)
                    nitrile_precursor = True
                    azide_precursor = True
                    # Add the co-occurrence constraint
                    if {"type": "co-occurrence", "details": {"targets": ["tetrazole", "Nitrile", "Azide"], "scope": "reaction_step", "description": "As an alternative to a named reaction, tetrazole formation is identified if a tetrazole ring is in the product, and both nitrile and azide functional groups are present in the reactants of the same reaction step."}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["tetrazole", "Nitrile", "Azide"], "scope": "reaction_step", "description": "As an alternative to a named reaction, tetrazole formation is identified if a tetrazole ring is in the product, and both nitrile and azide functional groups are present in the reactants of the same reaction step."}})

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

    # Check if tetrazole formation is late-stage (depth <= 2)
    is_late_stage = tetrazole_depth <= 2

    result = tetrazole_formed and nitrile_precursor and is_late_stage

    if result and is_late_stage:
        # Add the positional constraint if it's late-stage and the overall condition is met
        if {"type": "positional", "details": {"target": "ring_formation", "description": "The tetrazole formation event must occur at a synthetic depth of 2 or less.", "position_type": "depth", "operator": "<=", "value": 2}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ring_formation", "description": "The tetrazole formation event must occur at a synthetic depth of 2 or less.", "position_type": "depth", "operator": "<=", "value": 2}})

    return result, findings_json
