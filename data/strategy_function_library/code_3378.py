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


CROSS_COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki",
    "Negishi coupling",
    "Stille reaction_aryl",
    "Sonogashira alkyne_aryl halide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where a nitrile-containing component is used in a specific cross-coupling reaction. The function checks for the following reaction types: Suzuki coupling with boronic acids, Suzuki coupling with boronic esters, Suzuki, Negishi coupling, Stille reaction_aryl, and Sonogashira alkyne_aryl halide.
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

    # Initialize tracking variables
    has_aromatic_coupling = False
    has_nitrile_component = False
    has_biaryl_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_aromatic_coupling, has_nitrile_component, findings_json

        # Optimization: stop traversing if the strategy has already been found
        if has_aromatic_coupling:
            return

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            # Check if the reaction is a known cross-coupling type
            is_coupling_reaction = False
            detected_coupling_reaction = None
            for rxn_type in CROSS_COUPLING_REACTIONS:
                if checker.check_reaction(rxn_type, rsmi):
                    is_coupling_reaction = True
                    detected_coupling_reaction = rxn_type
                    break

            # If it's a coupling, check for a nitrile component in the same reaction
            if is_coupling_reaction:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                nitrile_is_present = False
                if any(checker.check_fg("Nitrile", r) for r in reactants) or \
                   checker.check_fg("Nitrile", product):
                    nitrile_is_present = True

                # If both conditions are met for this single reaction, set the flags and stop.
                if nitrile_is_present:
                    has_aromatic_coupling = True
                    has_nitrile_component = True

                    # Record findings
                    if detected_coupling_reaction and detected_coupling_reaction not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(detected_coupling_reaction)
                    if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

                    # Add structural constraint if both conditions are met in the same step
                    findings_json["structural_constraints"].append({
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "Nitrile",
                                [
                                    "Suzuki coupling with boronic acids",
                                    "Suzuki coupling with boronic esters",
                                    "Suzuki",
                                    "Negishi coupling",
                                    "Stille reaction_aryl",
                                    "Sonogashira alkyne_aryl halide"
                                ]
                            ],
                            "comment": "A Nitrile functional group must be present in the same reaction step as one of the specified cross-coupling reactions."
                        }
                    })
                    return

        # Continue traversing
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if both key elements were found in the same reaction step
    strategy_present = (has_aromatic_coupling or has_biaryl_formation) and has_nitrile_component
    return strategy_present, findings_json
