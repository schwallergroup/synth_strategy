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


TARGET_HETEROCYCLES = ["isoxazole", "oxadiazole", "isothiazole", "triazole"]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy where a nitrile group is converted
    to a specific heterocycle (isoxazole, oxadiazole, isothiazole, or triazole)
    via an oxime intermediate. The strategy is identified by finding a nitrile
    in an early stage of the synthesis, an oxime in a middle stage, and one of
    the target heterocycles in a late stage.
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

    # Track if we found the key functional groups and their depths
    nitrile_depth = None
    oxime_depth = None
    heterocycle_depth = None
    found_heterocycles = []

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_depth, oxime_depth, heterocycle_depth, findings_json, found_heterocycles

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for nitrile
            if checker.check_fg("Nitrile", mol_smiles):
                if nitrile_depth is None or depth > nitrile_depth:
                    nitrile_depth = depth
                    if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

            # Check for oxime (intermediate)
            if checker.check_fg("Oxime", mol_smiles):
                if oxime_depth is None or depth > oxime_depth:
                    oxime_depth = depth
                    if "Oxime" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Oxime")

            # Check for target heterocycles
            is_heterocycle_found_in_mol = False
            for h in TARGET_HETEROCYCLES:
                if checker.check_ring(h, mol_smiles):
                    is_heterocycle_found_in_mol = True
                    if heterocycle_depth is None or depth < heterocycle_depth:
                        heterocycle_depth = depth
                    if h not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(h)
                    if h not in found_heterocycles:
                        found_heterocycles.append(h)

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a 'mol' node
            next_depth = depth + 1

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    result = False
    # Check if we have the complete strategy: nitrile -> oxime -> heterocycle
    # This is identified by the chronological appearance of the functional groups
    # in the synthesis tree (high depth -> low depth).
    if (
        nitrile_depth is not None
        and oxime_depth is not None
        and heterocycle_depth is not None
        and nitrile_depth > oxime_depth > heterocycle_depth
    ):
        result = True
        # Add structural constraints if the full sequence is found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Nitrile",
                    "Oxime",
                    "target_heterocycle"
                ],
                "description": "The route must contain a Nitrile, an Oxime, and one of the target heterocycles (isoxazole, oxadiazole, isothiazole, or triazole)."
            }
        })
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_events": [
                    "Nitrile",
                    "Oxime",
                    "target_heterocycle"
                ],
                "description": "The entities must appear in a specific order based on synthetic depth, with Nitrile appearing earliest (highest depth), followed by Oxime, and finally the target heterocycle (lowest depth)."
            }
        })

    return result, findings_json
