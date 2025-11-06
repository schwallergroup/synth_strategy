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
    This function detects a synthetic strategy involving late-stage amide to ester conversion.
    Note: The route is traversed retrosynthetically, so we're looking for ester-to-amide
    conversions in the retrosynthetic direction.
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

    amide_to_ester_conversion = False
    max_depth_for_late_stage = 2  # Define what "late-stage" means in terms of depth

    def dfs_traverse(node, depth=0):
        nonlocal amide_to_ester_conversion, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Only consider late-stage reactions
            if depth <= max_depth_for_late_stage:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # In retrosynthesis, check for ester in reactants
                reactant_has_ester = False
                for r_smiles in reactants_smiles:
                    if checker.check_fg("Ester", r_smiles):
                        reactant_has_ester = True
                        if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ester")
                        break

                # Check for amide in product
                product_has_amide = False
                if checker.check_fg("Primary amide", product_smiles):
                    product_has_amide = True
                    if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                if checker.check_fg("Secondary amide", product_smiles):
                    product_has_amide = True
                    if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                if checker.check_fg("Tertiary amide", product_smiles):
                    product_has_amide = True
                    if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                if reactant_has_ester and product_has_amide:
                    # Add co-occurrence constraint if both are found
                    if {"type": "co-occurrence", "details": {"targets": ["Ester_in_retro_reactants", "Amide_in_retro_product"], "scope": "within_reaction"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Ester_in_retro_reactants", "Amide_in_retro_product"], "scope": "within_reaction"}})

                    # Additional check to ensure we're not just adding an amide group elsewhere
                    ester_count_reactants = sum(
                        1 for r in reactants_smiles if checker.check_fg("Ester", r)
                    )

                    ester_count_product = 1 if checker.check_fg("Ester", product_smiles) else 0

                    # In retrosynthesis, if ester count decreased and amide is present in product,
                    # it's likely a conversion from ester to amide (which represents amide to ester in forward synthesis)
                    if ester_count_reactants > ester_count_product:
                        amide_to_ester_conversion = True
                        # Add count constraint
                        if {"type": "count", "details": {"target": "Ester_in_retro_reactants", "operator": ">", "value_reference": "Ester_in_retro_product"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "Ester_in_retro_reactants", "operator": ">", "value_reference": "Ester_in_retro_product"}})
                        # Add positional constraint for late-stage
                        if {"type": "positional", "details": {"target": "amide_to_ester_conversion", "position": "late_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amide_to_ester_conversion", "position": "late_stage"}})
                        # Add named reaction if applicable
                        if "functional_group_interconversion" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("functional_group_interconversion")
                        return

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: amide_to_ester_conversion = {amide_to_ester_conversion}")
    return amide_to_ester_conversion, findings_json
