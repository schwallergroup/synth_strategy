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
    """This function detects the use of a Boc protection strategy. It identifies a multi-step sequence by checking for specific Boc protection and deprotection reactions, and the use of the resulting protected intermediate in subsequent steps. The specific named reactions checked for protection are defined in `BOC_PROTECTION_REACTIONS`, and for deprotection in `BOC_DEPROTECTION_REACTIONS`."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Track key elements of the protection strategy
    protected_molecules = {}  # Maps molecule SMILES to depth
    deprotected_molecules = {}  # Maps molecule SMILES to depth
    protection_reactions = []  # List of (depth, reactant, product) for protection reactions
    deprotection_reactions = []  # List of (depth, reactant, product) for deprotection reactions
    intermediate_reactions = []  # List of reactions using protected intermediates

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "reaction":
            if "rsmi" not in node["metadata"]:
                return

            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for Boc protection reaction
            protection_found = False
            for r in BOC_PROTECTION_REACTIONS:
                if checker.check_reaction(r, rsmi):
                    protection_found = True
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)
            
            if protection_found:
                print(f"Found Boc protection reaction at depth {depth}: {rsmi}")
                protection_reactions.append((depth, reactants_smiles, product_smiles))

                # Mark the product as Boc-protected
                if checker.check_fg("Boc", product_smiles):
                    if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Boc")
                    protected_molecules[product_smiles] = depth
                    print(f"Added protected molecule at depth {depth}: {product_smiles}")

            # Check for Boc deprotection reaction
            deprotection_found = False
            for r in BOC_DEPROTECTION_REACTIONS:
                if checker.check_reaction(r, rsmi):
                    deprotection_found = True
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)

            if deprotection_found:
                print(f"Found Boc deprotection reaction at depth {depth}: {rsmi}")
                deprotection_reactions.append((depth, reactants_smiles, product_smiles))

                # Mark reactants as being deprotected
                for reactant in reactants_smiles:
                    if checker.check_fg("Boc", reactant):
                        if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Boc")
                        deprotected_molecules[reactant] = depth
                        print(f"Added deprotected molecule at depth {depth}: {reactant}")

            # Check if any protected molecule is used in this reaction
            else:
                for reactant in reactants_smiles:
                    if checker.check_fg("Boc", reactant):
                        print(
                            f"Found reaction using protected intermediate at depth {depth}: {rsmi}"
                        )
                        intermediate_reactions.append((depth, reactants_smiles, product_smiles))
                        break

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Analyze the protection strategy
    has_protection = len(protection_reactions) > 0
    has_deprotection = len(deprotection_reactions) > 0
    has_intermediate_reactions = len(intermediate_reactions) > 0

    # Check if there's a proper sequence: protection -> intermediate reactions -> deprotection
    proper_sequence = False
    if has_protection and has_deprotection:
        # Get the minimum depth of protection and maximum depth of deprotection
        min_protection_depth = (
            min([depth for depth, _, _ in protection_reactions])
            if protection_reactions
            else float("inf")
        )
        max_deprotection_depth = (
            max([depth for depth, _, _ in deprotection_reactions]) if deprotection_reactions else 0
        )

        # In retrosynthetic analysis, higher depth means earlier in synthesis
        # So protection (earlier) should have higher depth than deprotection (later)
        if min_protection_depth > max_deprotection_depth:
            # Check if there are reactions between protection and deprotection
            intermediate_depths = [depth for depth, _, _ in intermediate_reactions]
            if (
                intermediate_depths
                and min(intermediate_depths) < min_protection_depth
                and max(intermediate_depths) > max_deprotection_depth
            ):
                proper_sequence = True
                print(
                    f"Found proper Boc protection sequence: protection at depth {min_protection_depth}, "
                    f"intermediate reactions, deprotection at depth {max_deprotection_depth}"
                )

    # A Boc protection strategy is present if we have either:
    # 1. A complete sequence (protection -> intermediate -> deprotection)
    # 2. Evidence of Boc-protected intermediates being used
    strategy_present = (
        has_protection and has_deprotection and proper_sequence
    ) or has_intermediate_reactions

    if strategy_present:
        print("Detected Boc protection strategy")
        # Add the structural constraint if the full strategy is detected
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_events": [
                    "Boc Protection",
                    "Intermediate Reaction with Boc Group",
                    "Boc Deprotection"
                ]
            }
        })
    else:
        if has_protection:
            print("Found Boc protection but not a complete protection strategy")
        if has_deprotection:
            print("Found Boc deprotection but not a complete protection strategy")

    return strategy_present, findings_json
