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

BOC_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a Boc protection/deprotection strategy by identifying reactions from the `BOC_PROTECTION_REACTIONS` and `BOC_DEPROTECTION_REACTIONS` lists. A strategy is flagged if at least one protection and one deprotection reaction are found, ideally with the deprotection occurring at a later synthetic stage (shallower depth) than the protection.
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

    # Initialize tracking variables
    protected_molecules = {}  # Track molecules that have been Boc protected
    deprotected_molecules = {}  # Track molecules that have been Boc deprotected

    # Check if any starting materials already contain Boc groups
    def check_starting_materials(node):
        if node["type"] == "mol" and node.get("in_stock", False):
            if checker.check_fg("Carbamic ester", node["smiles"]):
                print(f"Found starting material with Boc group: {node['smiles'][:20]}...")
                protected_molecules[node["smiles"]] = float(
                    "inf"
                )  # Mark as protected at "infinite" depth
                if "Carbamic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Carbamic ester")

        for child in node.get("children", []):
            check_starting_materials(child)

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "reaction":
            try:
                # Extract reaction SMILES and product SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                product_smiles = rsmi.split(">")[-1]

                # Check for Boc protection reaction using the predefined list
                for r_name in BOC_PROTECTION_REACTIONS:
                    if checker.check_reaction(r_name, rsmi):
                        print(f"Detected Boc protection reaction at depth {depth}")
                        protected_molecules[product_smiles] = depth
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        break # Found one, no need to check others for this reaction node

                # Check for Boc deprotection reaction using the predefined list
                for r_name in BOC_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(r_name, rsmi):
                        print(f"Detected Boc deprotection reaction at depth {depth}")
                        deprotected_molecules[product_smiles] = depth
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        break # Found one, no need to check others for this reaction node

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal to find starting materials with Boc groups
    check_starting_materials(route)

    # Start main traversal
    dfs_traverse(route)

    print(f"Protected molecules: {len(protected_molecules)}")
    print(f"Deprotected molecules: {len(deprotected_molecules)}")

    # Check if we have both protection and deprotection in the correct sequence
    has_sequence = False

    # Look for molecules that were protected and then deprotected
    for prot_mol, prot_depth in protected_molecules.items():
        for deprot_mol, deprot_depth in deprotected_molecules.items():
            # Deprotection should happen at a lower depth (later stage) than protection
            if deprot_depth < prot_depth:
                print(
                    f"Found protection at depth {prot_depth} and deprotection at depth {deprot_depth}"
                )
                has_sequence = True
                # Add structural constraint for sequence
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "before": "any_boc_protection_reaction",
                        "after": "any_boc_deprotection_reaction",
                        "description": "A Boc deprotection reaction must occur at a later synthetic stage (shallower depth) than a Boc protection reaction."
                    }
                })
                break
        if has_sequence:
            break

    # If we couldn't find a clear sequence but have both protection and deprotection,
    # still consider it a protection/deprotection strategy
    strategy_present = has_sequence or (
        len(protected_molecules) > 0 and len(deprotected_molecules) > 0
    )

    if not has_sequence and (len(protected_molecules) > 0 and len(deprotected_molecules) > 0):
        # Add structural constraint for co-occurrence if sequence not found but both present
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "any_boc_protection_reaction",
                    "any_boc_deprotection_reaction"
                ],
                "description": "The route must contain at least one Boc protection reaction and at least one Boc deprotection reaction. This serves as a fallback if a clear sequence is not identified."
            }
        })

    print(f"Protection/deprotection sequence detected: {strategy_present}")
    return strategy_present, findings_json