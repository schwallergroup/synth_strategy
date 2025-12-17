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


NITROGEN_FGS_OF_INTEREST = [
    "Oxime", "Primary amine", "Secondary amine", "Tertiary amine",
    "Primary amide", "Secondary amide", "Tertiary amide",
    "Imine", "Azide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a multi-part strategy where the synthetic route contains: (1) a ketone-to-secondary alcohol reduction step, (2) a secondary alcohol-to-ketone oxidation step, and (3) a ring-forming or ring-breaking step. Additionally, it verifies that the final product of the synthesis (at depth=1) contains both an aromatic halide and a nitrogen-containing functional group from a predefined list.
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

    # Track if we've found the key features
    found_ketone_to_alcohol = False
    found_alcohol_to_ketone = False
    found_ring_manipulation = False
    has_halogenated_aromatic = False
    ends_with_nitrogen_fg = False

    def dfs_traverse(node, depth):
        nonlocal found_ketone_to_alcohol, found_alcohol_to_ketone, found_ring_manipulation
        nonlocal has_halogenated_aromatic, ends_with_nitrogen_fg, findings_json

        if node.get("type") == "reaction":
            # Assume the reaction object is available on the node
            reaction = node["reaction"]

            # Check for ketone to alcohol transformation
            if checker.check_reactant_fg(reaction, "Ketone"):
                if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Ketone")
            if checker.check_product_fg(reaction, "Secondary alcohol"):
                if "Secondary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")

            if checker.check_reactant_fg(reaction, "Ketone") and \
               checker.check_product_fg(reaction, "Secondary alcohol"):
                found_ketone_to_alcohol = True
                if "ketone_to_alcohol_reduction" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ketone_to_alcohol_reduction")

            # Check for alcohol to ketone transformation
            if checker.check_reactant_fg(reaction, "Secondary alcohol"):
                if "Secondary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
            if checker.check_product_fg(reaction, "Ketone"):
                if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Ketone")

            if checker.check_reactant_fg(reaction, "Secondary alcohol") and \
               checker.check_product_fg(reaction, "Ketone"):
                found_alcohol_to_ketone = True
                if "alcohol_to_ketone_oxidation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("alcohol_to_ketone_oxidation")

            # Check for ring manipulation (any change in ring count)
            is_ring_forming = checker.is_ring_forming(reaction)
            is_ring_breaking = checker.is_ring_breaking(reaction)

            if is_ring_forming or is_ring_breaking:
                found_ring_manipulation = True
                if is_ring_forming and "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                if is_ring_breaking and "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")

            # In the final step (depth=1), check for specific features in the final product
            if depth == 1:
                if checker.check_product_fg(reaction, "Aromatic halide"):
                    has_halogenated_aromatic = True
                    if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

                nitrogen_fg_found_in_product = False
                for fg in NITROGEN_FGS_OF_INTEREST:
                    if checker.check_product_fg(reaction, fg):
                        nitrogen_fg_found_in_product = True
                        if fg not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg)
                if nitrogen_fg_found_in_product:
                    ends_with_nitrogen_fg = True

        # Traverse children, incrementing depth based on node type
        for child in node.get("children", []):
            new_depth = depth
            if node.get("type") != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root node at depth 1
    dfs_traverse(route, 1)

    strategy_present = (
        found_ketone_to_alcohol
        and found_alcohol_to_ketone
        and found_ring_manipulation
        and has_halogenated_aromatic
        and ends_with_nitrogen_fg
    )

    # Populate structural constraints based on the flags
    if found_ketone_to_alcohol:
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "ketone_to_alcohol_reduction", "operator": ">=", "value": 1, "description": "The route must contain at least one ketone-to-secondary alcohol reduction step. This is defined as a reaction where a 'Ketone' is a reactant and a 'Secondary alcohol' is a product."}})
    if found_alcohol_to_ketone:
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "alcohol_to_ketone_oxidation", "operator": ">=", "value": 1, "description": "The route must contain at least one secondary alcohol-to-ketone oxidation step. This is defined as a reaction where a 'Secondary alcohol' is a reactant and a 'Ketone' is a product."}})
    if found_ring_manipulation:
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": ["ring_formation", "ring_destruction"], "operator": ">=", "value": 1, "description": "The route must contain at least one step that is either ring-forming or ring-breaking."}})
    if has_halogenated_aromatic:
        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Aromatic halide", "position": "last_stage", "description": "The product of the final reaction step (depth=1) must contain an 'Aromatic halide' functional group."}})
    if ends_with_nitrogen_fg:
        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": ["Oxime", "Primary amine", "Secondary amine", "Tertiary amine", "Primary amide", "Secondary amide", "Tertiary amide", "Imine", "Azide"], "position": "last_stage", "description": "The product of the final reaction step (depth=1) must contain at least one of the specified nitrogen-containing functional groups."}})

    if strategy_present:
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["ketone_to_alcohol_reduction", "alcohol_to_ketone_oxidation", "ring_manipulation", "final_product_has_aromatic_halide", "final_product_has_nitrogen_fg"], "description": "The overall strategy requires the co-occurrence of five distinct features: a ketone reduction, an alcohol oxidation, a ring manipulation, and specific functional groups in the final product."}})

    return strategy_present, findings_json
