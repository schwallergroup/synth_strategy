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
    Detects a synthetic route that contains an isoxazole core and also includes at least one reaction step that appears to be an isoxazole halogenation or a cross-coupling from a halogenated isoxazole. This check does not enforce a specific sequence of these reactions.
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

    # Track key features
    has_isoxazole_halogenation = False
    has_isoxazole_coupling = False
    isoxazole_present = False

    def is_isoxazole_halogenation(reaction_node):
        """Check if a reaction involves halogenation of an isoxazole"""
        nonlocal findings_json
        if "rsmi" not in reaction_node["metadata"]:
            return False

        rsmi = reaction_node["metadata"]["mapped_reaction_smiles"]
        reactants = rsmi.split(">")[0].split(".")
        product = rsmi.split(">")[-1]

        # Check for isoxazole in reactants
        has_isoxazole_reactant = False
        for reactant in reactants:
            if reactant.strip():
                try:
                    if checker.check_ring("isoxazole", reactant):
                        has_isoxazole_reactant = True
                        if "isoxazole" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("isoxazole")
                        break
                except Exception:
                    continue

        # Check for halogenated isoxazole in product
        has_halogenated_product = False
        try:
            if checker.check_ring("isoxazole", product):
                if "isoxazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("isoxazole")
                # Check if product has halogen attached to isoxazole
                if checker.check_fg("Aromatic halide", product):
                    has_halogenated_product = True
                    if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                if checker.check_fg("Primary halide", product):
                    has_halogenated_product = True
                    if "Primary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary halide")
                if checker.check_fg("Secondary halide", product):
                    has_halogenated_product = True
                    if "Secondary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary halide")
                if checker.check_fg("Tertiary halide", product):
                    has_halogenated_product = True
                    if "Tertiary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary halide")
                if checker.check_fg("Alkenyl halide", product):
                    has_halogenated_product = True
                    if "Alkenyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Alkenyl halide")
        except Exception:
            pass

        if has_isoxazole_reactant and has_halogenated_product:
            if "isoxazole_halogenation" not in findings_json["atomic_checks"]["named_reactions"]:
                findings_json["atomic_checks"]["named_reactions"].append("isoxazole_halogenation")
        return has_isoxazole_reactant and has_halogenated_product

    def is_isoxazole_coupling(reaction_node):
        """Check if a reaction involves coupling with a halogenated isoxazole"""
        nonlocal findings_json
        if "rsmi" not in reaction_node["metadata"]:
            return False

        rsmi = reaction_node["metadata"]["mapped_reaction_smiles"]
        reactants = rsmi.split(">")[0].split(".")
        product = rsmi.split(">")[-1]

        # Check for halogenated isoxazole in reactants
        has_halogenated_isoxazole_reactant = False
        for reactant in reactants:
            if reactant.strip():
                try:
                    if checker.check_ring("isoxazole", reactant):
                        if "isoxazole" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("isoxazole")
                        if checker.check_fg("Aromatic halide", reactant):
                            has_halogenated_isoxazole_reactant = True
                            if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                        if checker.check_fg("Primary halide", reactant):
                            has_halogenated_isoxazole_reactant = True
                            if "Primary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary halide")
                        if checker.check_fg("Secondary halide", reactant):
                            has_halogenated_isoxazole_reactant = True
                            if "Secondary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary halide")
                        if checker.check_fg("Tertiary halide", reactant):
                            has_halogenated_isoxazole_reactant = True
                            if "Tertiary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Tertiary halide")
                        if checker.check_fg("Alkenyl halide", reactant):
                            has_halogenated_isoxazole_reactant = True
                            if "Alkenyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Alkenyl halide")
                        if has_halogenated_isoxazole_reactant:
                            break
                except Exception:
                    continue

        # Check for coupled isoxazole in product
        has_coupled_product = False
        has_isoxazole_product = False
        try:
            if checker.check_ring("isoxazole", product):
                has_isoxazole_product = True
                if "isoxazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("isoxazole")
                # For coupling reactions, we expect the halogen to be replaced
                # with a carbon-based group, so we check if the halogen is gone
                if not (
                    checker.check_fg("Aromatic halide", product)
                    or checker.check_fg("Primary halide", product)
                    or checker.check_fg("Secondary halide", product)
                    or checker.check_fg("Tertiary halide", product)
                    or checker.check_fg("Alkenyl halide", product)
                ):
                    has_coupled_product = True
        except Exception:
            pass

        if has_halogenated_isoxazole_reactant and has_isoxazole_product and has_coupled_product:
            if "isoxazole_cross_coupling" not in findings_json["atomic_checks"]["named_reactions"]:
                findings_json["atomic_checks"]["named_reactions"].append("isoxazole_cross_coupling")
        return has_halogenated_isoxazole_reactant and has_isoxazole_product and has_coupled_product

    def has_isoxazole(smiles):
        """Check if molecule contains an isoxazole core"""
        nonlocal findings_json
        try:
            result = checker.check_ring("isoxazole", smiles)
            if result:
                if "isoxazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("isoxazole")
            return result
        except Exception:
            return False

    def dfs_traverse(node, depth=0):
        nonlocal has_isoxazole_halogenation, has_isoxazole_coupling, isoxazole_present, findings_json

        # Check for isoxazole in molecules
        if node["type"] == "mol" and "smiles" in node:
            if has_isoxazole(node["smiles"]):
                isoxazole_present = True

        # Check for isoxazole halogenation and coupling
        if node["type"] == "reaction" and "metadata" in node:
            if is_isoxazole_halogenation(node):
                has_isoxazole_halogenation = True

            if is_isoxazole_coupling(node):
                has_isoxazole_coupling = True

        # Recursively process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Determine if this matches our target strategy
    strategy_detected = isoxazole_present and (has_isoxazole_halogenation or has_isoxazole_coupling)

    # Populate structural constraints based on detected flags
    if isoxazole_present:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "isoxazole",
                "operator": ">=",
                "value": 1
            }
        })
    if has_isoxazole_halogenation or has_isoxazole_coupling:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target_group": [
                    "isoxazole_halogenation",
                    "isoxazole_cross_coupling"
                ],
                "operator": ">=",
                "value": 1
            }
        })

    return strategy_detected, findings_json
