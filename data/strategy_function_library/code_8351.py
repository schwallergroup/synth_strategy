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


# Refactoring for Enumeration: Isolate the lists of reaction names
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
    This function detects a Boc protection/deprotection sequence in the synthesis.
    Returns True if a Boc group is added early in the synthesis and removed later.
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

    # Track protection and deprotection events with depth information
    protection_events = []  # Will store (depth, molecule_smiles, atom_indices)
    deprotection_events = []  # Will store (depth, molecule_smiles, atom_indices)

    def dfs_traverse(node, depth=0):
        nonlocal findings_json, protection_events, deprotection_events
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc protection using checker function
            for r_name in BOC_PROTECTION_REACTIONS:
                if checker.check_reaction(r_name, rsmi):
                    print(f"Found Boc protection reaction at depth {depth}: {rsmi}")
                    # Store the protection event with depth and product molecule
                    protection_events.append((depth, product))
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    break # Found a protection reaction, no need to check others

            # Check for Boc deprotection using checker function
            for r_name in BOC_DEPROTECTION_REACTIONS:
                if checker.check_reaction(r_name, rsmi):
                    print(f"Found Boc deprotection reaction at depth {depth}: {rsmi}")
                    # Store the deprotection event with depth and reactant molecule
                    # In deprotection, the Boc group is on the reactant
                    for reactant in reactants:
                        if checker.check_fg("Carbamic ester", reactant):
                            deprotection_events.append((depth, reactant))
                            if "Carbamic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Carbamic ester")
                            break # Found the FG in reactant, no need to check other reactants
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    break # Found a deprotection reaction, no need to check others

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is a chemical node
            next_depth = depth + 1

        # Traverse children (going backward in synthesis)
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if we have both protection and deprotection events
    has_protection = len(protection_events) > 0
    has_deprotection = len(deprotection_events) > 0

    print(f"Boc protection events: {len(protection_events)}")
    print(f"Boc deprotection events: {len(deprotection_events)}")

    # Check if protection happens at a higher depth (earlier in synthesis)
    # than deprotection (later in synthesis)
    valid_sequence = False
    if has_protection and has_deprotection:
        # Get the minimum depth of protection (earliest protection)
        min_protection_depth = min(depth for depth, _ in protection_events)
        # Get the maximum depth of deprotection (latest deprotection)
        max_deprotection_depth = max(depth for depth, _ in deprotection_events)

        # In a valid sequence, protection should happen earlier (higher depth)
        # than deprotection (lower depth)
        valid_sequence = min_protection_depth > max_deprotection_depth
        print(f"Min protection depth: {min_protection_depth}")
        print(f"Max deprotection depth: {max_deprotection_depth}")
        print(f"Valid sequence: {valid_sequence}")

    result = has_protection and has_deprotection and valid_sequence

    if has_protection and has_deprotection:
        # Co-occurrence constraint met
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Boc_protection",
                    "Boc_deprotection"
                ],
                "description": "The route must contain at least one Boc protection reaction and at least one Boc deprotection reaction."
            }
        })

    if valid_sequence:
        # Sequence constraint met
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_events": [
                    "Boc_protection",
                    "Boc_deprotection"
                ],
                "description": "A Boc protection reaction must occur earlier in the synthesis (higher depth) than a Boc deprotection reaction."
            }
        })

    # Return True if we have both protection and deprotection in the correct sequence
    return result, findings_json
