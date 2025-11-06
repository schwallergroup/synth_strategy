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
    This function detects a silyl protection-deprotection sequence in the synthesis.
    Specifically looking for O-Si bond formation followed by cleavage.
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

    # Track protection and deprotection events
    protection_events = []  # Will store (depth, product_smiles, reactant_smiles)
    deprotection_events = []  # Will store (depth, product_smiles, reactant_smiles)

    result = False # Initialize the main boolean result

    def dfs_traverse(node, depth=0):
        nonlocal protection_events, deprotection_events, findings_json
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Extract reactants and product
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]
                print(f"  Product: {product}")
                print(f"  Reactants: {reactants}")
            except Exception as e:
                print(f"Error extracting reactants/product: {e}")
                return

            # Check for silyl protection reaction - either by reaction type or by functional group changes
            if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                print(f"Silyl protection found at depth {depth} (by reaction type)")
                protection_events.append((depth, product, reactants))
                if "Alcohol protection with silyl ethers" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Alcohol protection with silyl ethers")
            # Alternative detection: check if alcohol in reactants and silyl group in product
            elif any(checker.check_fg("Alcohol", r) for r in reactants) and checker.check_fg(
                "Silyl protective group", product
            ):
                # Verify this is likely a protection by checking for silyl chloride in reactants
                if any("Si" in r and "Cl" in r for r in reactants):
                    print(f"Silyl protection found at depth {depth} (by functional group analysis)")
                    protection_events.append((depth, product, reactants))
                    if "Alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Alcohol")
                    if "Silyl protective group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Silyl protective group")

            # Check for silyl deprotection reactions - either by reaction type or by functional group changes
            deprotection_reaction_found = False
            if (
                checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
            ):
                print(f"Silyl deprotection found at depth {depth} (by reaction type)")
                deprotection_events.append((depth, product, reactants))
                if "Alcohol deprotection from silyl ethers" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Alcohol deprotection from silyl ethers")
                deprotection_reaction_found = True
            if (
                checker.check_reaction("Alcohol deprotection from silyl ethers (double)", rsmi)
            ):
                print(f"Silyl deprotection found at depth {depth} (by reaction type)")
                deprotection_events.append((depth, product, reactants))
                if "Alcohol deprotection from silyl ethers (double)" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Alcohol deprotection from silyl ethers (double)")
                deprotection_reaction_found = True
            if (
                checker.check_reaction("Alcohol deprotection from silyl ethers (diol)", rsmi)
            ):
                print(f"Silyl deprotection found at depth {depth} (by reaction type)")
                deprotection_events.append((depth, product, reactants))
                if "Alcohol deprotection from silyl ethers (diol)" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Alcohol deprotection from silyl ethers (diol)")
                deprotection_reaction_found = True

            # Alternative detection: check if silyl group in reactants and alcohol in product
            if not deprotection_reaction_found and any(
                checker.check_fg("Silyl protective group", r) for r in reactants
            ) and checker.check_fg("Alcohol", product):
                print(f"Silyl deprotection found at depth {depth} (by functional group analysis)")
                deprotection_events.append((depth, product, reactants))
                if "Silyl protective group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Silyl protective group")
                if "Alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Alcohol")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    print(f"Protection events: {len(protection_events)}")
    print(f"Deprotection events: {len(deprotection_events)}")

    # Check if both protection and deprotection were found
    if not protection_events or not deprotection_events:
        print("Missing either protection or deprotection events")
        return False, findings_json

    # In retrosynthesis, higher depth means earlier stage in forward synthesis
    # Protection should occur before deprotection in forward synthesis
    for prot_depth, prot_product, prot_reactants in protection_events:
        for deprot_depth, deprot_product, deprot_reactants in deprotection_events:
            # Check if protection occurs before deprotection in forward synthesis
            if prot_depth > deprot_depth:
                print(
                    f"Found potential sequence: protection at depth {prot_depth}, deprotection at depth {deprot_depth}"
                )

                # Check if the protected molecule is related to the deprotected one
                # In a proper sequence, the product of protection should be involved in deprotection
                for reactant in deprot_reactants:
                    if reactant == prot_product:
                        print("Silyl protection-deprotection strategy detected")
                        result = True
                        # Add structural constraints if the sequence is detected
                        if {"type": "co-occurrence", "details": {"targets": ["silyl_protection_event", "silyl_deprotection_event"]}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["silyl_protection_event", "silyl_deprotection_event"]}})
                        if {"type": "sequence", "details": {"ordered_events": ["silyl_protection_event", "silyl_deprotection_event"], "molecule_flow_constraint": "product_of_first_is_reactant_of_second"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["silyl_protection_event", "silyl_deprotection_event"], "molecule_flow_constraint": "product_of_first_is_reactant_of_second"}})
                        return result, findings_json

    print("No valid silyl protection-deprotection sequence found")
    return result, findings_json
