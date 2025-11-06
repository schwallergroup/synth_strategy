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


BOC_PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a strategy where Boc protection of a nitrogen
    (typically in a piperidine) occurs as the final step of the synthesis.
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

    found_boc_protection = False
    min_boc_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal found_boc_protection, min_boc_depth, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants_smiles = parts[0].split(".")
            product_smiles = parts[2]

            # Check if this is a Boc protection reaction using a predefined list
            is_boc_protection = False
            for name in BOC_PROTECTION_REACTIONS:
                if checker.check_reaction(name, rsmi):
                    is_boc_protection = True
                    if name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                    break

            # If reaction check fails, try to determine by functional group changes
            if not is_boc_protection:
                # Check if product contains Boc group
                has_boc_in_product = checker.check_fg("Boc", product_smiles)
                if has_boc_in_product and "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boc")

                if has_boc_in_product:
                    # Check if any reactant has an unprotected nitrogen
                    for reactant in reactants_smiles:
                        has_boc_in_reactant = checker.check_fg("Boc", reactant)

                        # Check for various types of amines that can be protected
                        has_amine = False
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

                        # Check for piperidine specifically
                        has_piperidine = checker.check_ring("piperidine", reactant)
                        if has_piperidine and "piperidine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("piperidine")

                        if (has_amine or has_piperidine) and not has_boc_in_reactant:
                            is_boc_protection = True
                            break

            if is_boc_protection:
                found_boc_protection = True
                min_boc_depth = min(min_boc_depth, depth)

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)

    # Consider it late-stage if it occurs at depth 0, 1, 2, or 3
    result = found_boc_protection and min_boc_depth <= 3

    if result:
        # Add the structural constraint if the overall condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Boc protection",
                "position_type": "min_depth",
                "operator": "<=",
                "value": 3
            }
        })

    return result, findings_json
