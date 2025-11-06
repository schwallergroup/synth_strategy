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


HALIDE_CLASSES_TO_CHECK = [
    "Aromatic halide",
    "Primary halide",
    "Secondary halide",
    "Tertiary halide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects whether CF3 and chloro substituents are retained throughout
    the synthesis, indicating a strategy that preserves these functional groups.
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

    # First, find the final product (root of the tree)
    final_product = route if route["type"] == "mol" else None

    if not final_product:
        return False, findings_json

    def _is_chloro_present(smiles_string):
        if "Cl" not in smiles_string:
            return False
        found_chloro = False
        for hc in HALIDE_CLASSES_TO_CHECK:
            if checker.check_fg(hc, smiles_string):
                if hc not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append(hc)
                found_chloro = True
        return found_chloro

    # Check if final product has CF3 and chloro groups
    has_cf3 = checker.check_fg("Trifluoro group", final_product["smiles"])
    if has_cf3:
        if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
            findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

    has_chloro = _is_chloro_present(final_product["smiles"])

    result = True

    if not (has_cf3 and has_chloro):
        result = False
    else:
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Trifluoro group", "chloro_group"], "scope": "final_product"}})

    # Track if these groups are introduced during synthesis
    cf3_introduced = False
    chloro_introduced = False

    def check_reactions(node, depth=0):
        nonlocal cf3_introduced, chloro_introduced, findings_json

        # Check reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if reactants have CF3 and chloro
            reactants_have_cf3 = any(
                checker.check_fg("Trifluoro group", r) for r in reactants_smiles
            )
            if reactants_have_cf3:
                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

            reactants_have_chloro = any(_is_chloro_present(r) for r in reactants_smiles)

            # Check if product has CF3 and chloro
            product_has_cf3 = checker.check_fg("Trifluoro group", product_smiles)
            if product_has_cf3:
                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

            product_has_chloro = _is_chloro_present(product_smiles)

            # Check if CF3 or chloro is introduced in this reaction
            if product_has_cf3 and not reactants_have_cf3:
                cf3_introduced = True

            if product_has_chloro and not reactants_have_chloro:
                chloro_introduced = True

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            check_reactions(child, new_depth)

    # Start traversal from the root
    check_reactions(route)

    # Strategy is present if both groups are in final product and NOT introduced during synthesis
    if cf3_introduced:
        result = False
    else:
        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "introduction_of_Trifluoro group", "scope": "synthesis_route"}})

    if chloro_introduced:
        result = False
    else:
        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "introduction_of_chloro_group", "scope": "synthesis_route"}})

    return result, findings_json
