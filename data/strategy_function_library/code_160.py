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
    This function detects if the synthesis involves a late-stage esterification
    with a complex fragment containing a trifluoromethoxy group.
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

    esterification_detected = False
    esterification_found_at_late_stage = False
    trifluoromethoxy_found_in_reactant = False
    complex_fragment_found_with_trifluoromethoxy = False
    complex_fragment_found_in_reactants = False

    def is_complex_fragment(smiles):
        """Check if a molecule is a complex fragment (has multiple rings or functional groups)"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Check for complexity: rings
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() >= 1:
            print(f"Complex due to rings: {ring_info.NumRings()} rings")
            # No specific ring system name to add, but indicates complexity
            return True

        # Count different functional groups
        fg_count = 0
        for fg in [
            "Aromatic halide",
            "Ether",
            "Ester",
            "Amide",
            "Nitrile",
            "Nitro group",
            "Carboxylic acid",
            "Phenol",
            "Alcohol",
            "Trifluoro group",
        ]:
            if checker.check_fg(fg, smiles):
                fg_count += 1
                if fg not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append(fg)
                print(f"Found functional group: {fg}")

        if fg_count >= 2:  # At least 2 functional groups for complexity
            print(f"Complex due to functional groups: {fg_count}")
            if {"type": "count", "details": {"target": "functional_groups_on_complex_fragment", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "count", "details": {"target": "functional_groups_on_complex_fragment", "operator": ">=", "value": 2}})
            return True

        # If molecule has at least 12 atoms, consider it complex
        if mol.GetNumHeavyAtoms() >= 12:
            print(f"Complex due to size: {mol.GetNumHeavyAtoms()} heavy atoms")
            return True

        return False

    def has_trifluoromethoxy(smiles):
        """Check if a molecule contains a trifluoromethoxy group"""
        nonlocal trifluoromethoxy_found_in_reactant
        # Use checker to detect trifluoro group
        if checker.check_fg("Trifluoro group", smiles):
            if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
            print(f"Found trifluoro group in: {smiles}")

            # Check if it's connected to oxygen (trifluoromethoxy)
            if "OC(F)(F)F" in smiles or "OCF3" in smiles:
                print(f"Found trifluoromethoxy pattern in SMILES: {smiles}")
                trifluoromethoxy_found_in_reactant = True
                return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal esterification_detected, esterification_found_at_late_stage, complex_fragment_found_with_trifluoromethoxy, complex_fragment_found_in_reactants

        if node["type"] == "reaction":
            if depth <= 2:  # Late stage (low depth)
                if {"type": "positional", "details": {"target": "Esterification", "position": "late_stage", "max_depth": 2}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Esterification", "position": "late_stage", "max_depth": 2}})
                esterification_found_at_late_stage = True

            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for esterification reactions - try multiple types
                is_esterification = False
                esterification_reaction_types = [
                    "Esterification of Carboxylic Acids",
                    "Transesterification",
                    "Schotten-Baumann to ester",
                    "O-alkylation of carboxylic acids with diazo compounds",
                    "Oxidative esterification of primary alcohols",
                ]
                for rxn_type in esterification_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Found esterification reaction type: {rxn_type}")
                        is_esterification = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if is_esterification:
                    # Extract reactants
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if any reactant has trifluoromethoxy group and is complex
                    for r in reactants:
                        if has_trifluoromethoxy(r):
                            print(f"Found reactant with trifluoromethoxy group: {r}")
                            if is_complex_fragment(r):
                                print(
                                    f"Confirmed complex fragment with trifluoromethoxy group: {r}"
                                )
                                complex_fragment_found_with_trifluoromethoxy = True
                                esterification_detected = True
                                # No early return here to allow all findings to be collected

                    # Also check if the product contains a trifluoromethoxy group
                    # This handles cases where the trifluoromethoxy group is in the product
                    if has_trifluoromethoxy(product):
                        print(f"Found product with trifluoromethoxy group: {product}")
                        # Check if any reactant is complex
                        for r in reactants:
                            if is_complex_fragment(r):
                                print(f"Confirmed complex fragment in reactants: {r}")
                                complex_fragment_found_in_reactants = True
                                esterification_detected = True
                                # No early return here to allow all findings to be collected

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node["type"] != "reaction": # Current node is 'chemical'
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Final check for structural constraints based on collected flags
    if esterification_found_at_late_stage and (complex_fragment_found_with_trifluoromethoxy or complex_fragment_found_in_reactants):
        if {"type": "co-occurrence", "details": {"targets": ["Esterification", "Trifluoro group"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Esterification", "Trifluoro group"]}})

    print(f"Late-stage esterification with complex fragment: {esterification_detected}")
    return esterification_detected, findings_json
