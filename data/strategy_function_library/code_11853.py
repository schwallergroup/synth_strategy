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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy where a nitrile group is preserved throughout
    the synthesis and a tertiary alcohol is formed via addition to a ketone.
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

    nitrile_present = False
    tertiary_alcohol_formed = False
    ketone_present = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_present, tertiary_alcohol_formed, ketone_present, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check if this reaction forms a tertiary alcohol from a ketone
            if checker.check_fg("Ketone", reactants):
                if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Ketone")
            if checker.check_fg("Tertiary alcohol", product):
                if "Tertiary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")

            if checker.check_fg("Ketone", reactants) and checker.check_fg(
                "Tertiary alcohol", product
            ):
                tertiary_alcohol_formed = True
                ketone_present = True
                print(
                    f"Depth {depth}: Tertiary alcohol formation from ketone detected in reaction: {rsmi}"
                )

            # Verify nitrile preservation across the reaction
            nitrile_in_reactants = checker.check_fg("Nitrile", reactants)
            nitrile_in_product = checker.check_fg("Nitrile", product)

            if nitrile_in_reactants:
                if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
            if nitrile_in_product:
                if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

            if nitrile_in_reactants and nitrile_in_product:
                nitrile_present = True
                print(f"Depth {depth}: Nitrile preserved in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if nitrile is preserved and tertiary alcohol is formed from ketone
    result = nitrile_present and tertiary_alcohol_formed and ketone_present
    if result:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "nitrile_preservation",
                    "tertiary_alcohol_formation_from_ketone"
                ],
                "description": "The overall strategy requires both the preservation of a nitrile group in at least one reaction step and the formation of a tertiary alcohol from a ketone in at least one reaction step."
            }
        })

    print(
        f"Nitrile preserved: {nitrile_present}, Tertiary alcohol formed: {tertiary_alcohol_formed}, Ketone present: {ketone_present}"
    )
    return result, findings_json
