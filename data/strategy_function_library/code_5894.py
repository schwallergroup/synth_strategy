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
    This function detects a synthetic strategy where a primary alcohol is oxidized to an ester
    in the final step while maintaining nitrogen protection throughout the synthesis.
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

    # Track if we found the pattern
    found_alcohol_oxidation = False
    found_boc_protection = False

    # Track functional groups at each step
    molecules_with_boc = []

    def dfs_traverse(node, depth=0):
        nonlocal found_alcohol_oxidation, found_boc_protection, findings_json

        if node["type"] == "mol":
            # Check if molecule has BOC protection
            if checker.check_fg("Boc", node["smiles"]):
                molecules_with_boc.append({"depth": depth, "smiles": node["smiles"]})
                if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boc")
                # print(f"Found Boc protection at depth {depth}: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a reaction involving alcohol oxidation
            # Find reactant with primary alcohol
            reactant_with_alcohol = None
            for r in reactants:
                if checker.check_fg("Primary alcohol", r):
                    reactant_with_alcohol = r
                    if "Primary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                    break

            # Check if product has ester
            has_ester = checker.check_fg("Ester", product)
            if has_ester:
                if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Ester")

            if reactant_with_alcohol and has_ester:
                # Check for various oxidation reactions
                oxidation_reactions = [
                    "Oxidative esterification of primary alcohols",
                    "Oxidation of alcohol and aldehyde to ester",
                ]

                for rxn_type in oxidation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        found_alcohol_oxidation = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        # print(f"Found alcohol oxidation ({rxn_type}) at depth {depth}: {rsmi}")
                        break

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if we have Boc protection throughout the synthesis
    if molecules_with_boc:
        # Get all depth levels
        depths = set(mol["depth"] for mol in molecules_with_boc)

        if depths:
            min_depth = min(depths)
            max_depth = max(depths)

            # Check if Boc protection is maintained throughout the synthesis
            # We need Boc protection at the final product (depth 0) and at least one earlier step
            if 0 in depths:
                # Structural constraint: Boc at last stage
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Boc", "position": "last_stage"}})
                if max_depth > 0:
                    found_boc_protection = True
                    # Structural constraint: Boc at not last stage
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Boc", "position": "not_last_stage"}})
                    # print(f"Found Boc protection at depths: {sorted(depths)}")

    # Check if we have the complete pattern
    result = found_alcohol_oxidation and found_boc_protection
    if result:
        # Structural constraint: co-occurrence of Boc, Primary alcohol, Ester
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Boc", "Primary alcohol", "Ester"]}})
        # print("Strategy detected: Late-stage alcohol oxidation with preserved N-protection")
    # else:
        # print(
            # f"Strategy not detected: alcohol_oxidation={found_alcohol_oxidation}, boc_protection={found_boc_protection}"
        # )
    return result, findings_json