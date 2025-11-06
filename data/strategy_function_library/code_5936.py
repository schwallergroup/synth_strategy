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


LATE_STAGE_ETHERIFICATIONS = [
    "Williamson Ether Synthesis",
    "Williamson Ether Synthesis (intra to epoxy)",
    "Mitsunobu aryl ether",
    "Alcohol to ether",
    "Chan-Lam etherification",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a build-connect-diversify strategy where core fragments are connected
    through epoxide chemistry, followed by late-stage diversification through etherification.
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

    # We'll track the presence of each component of the strategy
    has_epoxide_connection = False
    has_late_diversification = False
    reaction_sequence = []

    def dfs_traverse(node, depth=0):
        nonlocal has_epoxide_connection, has_late_diversification, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Skip if we can't process the reaction
            if not reactants or not product:
                return

            # Check for epoxide connection (typically in earlier steps)
            # Look for epoxide ring-opening reactions
            epoxide_reactant = False
            for reactant in reactants:
                if checker.check_ring("oxirane", reactant):
                    epoxide_reactant = True
                    findings_json["atomic_checks"]["ring_systems"].append("oxirane")
                    break

            # Verify the epoxide is being opened (not present in product)
            if epoxide_reactant and not checker.check_ring("oxirane", product):
                has_epoxide_connection = True
                reaction_sequence.append((depth, "epoxide_connection"))
                findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")
                print(f"Found epoxide connection at depth {depth}")

            # Check for late-stage diversification through etherification
            # Allow for depth 0, 1, or 2 to be considered late-stage
            if depth <= 2:
                ether_formed_by_fg = checker.check_fg("Ether", product) and not all(
                    checker.check_fg("Ether", r) for r in reactants
                )
                if ether_formed_by_fg:
                    findings_json["atomic_checks"]["functional_groups"].append("Ether")
                    findings_json["atomic_checks"]["named_reactions"].append("functional_group_formation")

                ether_formed_by_name = False
                for name in LATE_STAGE_ETHERIFICATIONS:
                    if checker.check_reaction(name, rsmi):
                        ether_formed_by_name = True
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                        break

                if ether_formed_by_fg or ether_formed_by_name:
                    has_late_diversification = True
                    reaction_sequence.append((depth, "etherification"))
                    # Add positional constraint if met
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "etherification",
                            "position": "late_stage (depth <= 2)"
                        }
                    })
                    print(
                        f"Found late-stage diversification through etherification at depth {depth}"
                    )

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the sequence is in the correct order (epoxide connection before etherification)
    correct_sequence = True
    if reaction_sequence:
        epoxide_depths = [d for d, r_type in reaction_sequence if r_type == "epoxide_connection"]
        ether_depths = [d for d, r_type in reaction_sequence if r_type == "etherification"]
        if epoxide_depths and ether_depths:
            # Epoxide connection should be at greater depth (earlier in synthesis)
            if min(epoxide_depths) <= max(ether_depths):
                correct_sequence = False
                print("Incorrect sequence: epoxide connection should occur before etherification")
            else:
                # Add sequence constraint if met
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "before": "epoxide_ring_opening",
                        "after": "etherification"
                    }
                })

    # Check if we have the complete strategy
    result = False
    if has_epoxide_connection and has_late_diversification and correct_sequence:
        result = True
        # Add co-occurrence constraint if met
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "epoxide_ring_opening",
                    "late_stage_etherification"
                ]
            }
        })
        print(
            "Found build-connect-diversify strategy with epoxide connection and late-stage etherification"
        )
    else:
        print(
            f"Strategy incomplete: epoxide_connection={has_epoxide_connection}, late_diversification={has_late_diversification}, correct_sequence={correct_sequence}"
        )

    return result, findings_json
