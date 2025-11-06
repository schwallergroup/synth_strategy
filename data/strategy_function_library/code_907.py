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


def main(route) -> Tuple[bool, Dict]:
    """
    Traverses a synthesis tree to detect a nitrobenzenesulfonyl (e.g., nosyl) protecting group strategy for amines. The strategy is flagged as present if the route contains at least one reaction for the protection (e.g., amine + nitrobenzenesulfonyl chloride -> nitrobenzenesulfonamide) and at least one reaction for the deprotection (nitrobenzenesulfonamide -> amine).
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

    nosyl_protected = False
    nosyl_deprotected = False

    def dfs_traverse(node, depth=0):
        nonlocal nosyl_protected, nosyl_deprotected, findings_json

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # --- Check for Nitrobenzenesulfonyl Deprotection ---
                # A reactant is a nitro-sulfonamide, and the product is an amine without the sulfonamide.
                reactant_is_protected = False
                for r in reactants:
                    if checker.check_fg("Sulfonamide", r):
                        if "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")
                    if checker.check_fg("Nitro group", r):
                        if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

                    if checker.check_fg("Sulfonamide", r) and checker.check_fg("Nitro group", r):
                        reactant_is_protected = True
                        break
                
                if reactant_is_protected:
                    product_is_amine = (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                    )
                    if checker.check_fg("Primary amine", product):
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    if checker.check_fg("Secondary amine", product):
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                    if checker.check_fg("Tertiary amine", product):
                        if "Tertiary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")

                    if product_is_amine and not checker.check_fg("Sulfonamide", product):
                        nosyl_deprotected = True
                        if "nosyl_deprotection" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("nosyl_deprotection")

                # --- Check for Nitrobenzenesulfonyl Protection ---
                # Reactants are an amine and a nitro-sulfonyl halide, product is a nitro-sulfonamide.
                has_amine_reactant = False
                has_sulfonyl_halide_reactant = False
                for r in reactants:
                    if checker.check_fg("Primary amine", r):
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        has_amine_reactant = True
                    if checker.check_fg("Secondary amine", r):
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                        has_amine_reactant = True
                    if checker.check_fg("Sulfonyl halide", r):
                        if "Sulfonyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Sulfonyl halide")
                    if checker.check_fg("Nitro group", r):
                        if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

                    if checker.check_fg("Sulfonyl halide", r) and checker.check_fg("Nitro group", r):
                        has_sulfonyl_halide_reactant = True
                
                if has_amine_reactant and has_sulfonyl_halide_reactant:
                    if checker.check_fg("Sulfonamide", product):
                        if "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")
                    if checker.check_fg("Nitro group", product):
                        if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

                    if checker.check_fg("Sulfonamide", product) and checker.check_fg("Nitro group", product):
                        nosyl_protected = True
                        if "nosyl_protection" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("nosyl_protection")

            except Exception:
                # Ignore reactions that cause parsing errors, as checkers might fail on invalid SMILES.
                pass

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., chemical), increase depth
                new_depth = depth + 1
            # If current node is a reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    if route:
        dfs_traverse(route)

    result = nosyl_protected and nosyl_deprotected

    if result:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "nosyl_protection",
                    "nosyl_deprotection"
                ]
            }
        })

    return result, findings_json
