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
    Detects a synthetic strategy involving sequential heteroatom alkylations (O then N)
    with a halogen exchange activation step between the alkylations.
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

    # Track if we've found each key step
    o_alkylation_found = False
    halogen_exchange_found = False
    n_alkylation_found = False

    # Track the order of reactions
    reaction_sequence = []

    def identify_reaction_type(rsmi):
        nonlocal findings_json
        # Extract reactants and products
        reactants_str = rsmi.split(">")[0]
        products_str = rsmi.split(">")[-1]
        reactant_mols = reactants_str.split(".")

        # Check for O-alkylation (Williamson Ether Synthesis or similar)
        if (
            checker.check_reaction("Williamson Ether Synthesis", rsmi)
            or checker.check_reaction("Williamson Ether Synthesis (intra to epoxy)", rsmi)
            or checker.check_reaction("Mitsunobu aryl ether", rsmi)
        ):
            if "Williamson Ether Synthesis" not in findings_json["atomic_checks"]["named_reactions"]:
                findings_json["atomic_checks"]["named_reactions"].append("Williamson Ether Synthesis")
            if "Williamson Ether Synthesis (intra to epoxy)" not in findings_json["atomic_checks"]["named_reactions"]:
                findings_json["atomic_checks"]["named_reactions"].append("Williamson Ether Synthesis (intra to epoxy)")
            if "Mitsunobu aryl ether" not in findings_json["atomic_checks"]["named_reactions"]:
                findings_json["atomic_checks"]["named_reactions"].append("Mitsunobu aryl ether")
            print("O-alkylation identified: Named reaction check")
            return "O-alkylation"

        # Check for halogen exchange (Finkelstein or similar)
        if checker.check_reaction("Finkelstein reaction", rsmi):
            if "Finkelstein reaction" not in findings_json["atomic_checks"]["named_reactions"]:
                findings_json["atomic_checks"]["named_reactions"].append("Finkelstein reaction")
            print("Halogen exchange identified: Finkelstein reaction")
            return "halogen-exchange"

        # Check for N-alkylation
        if checker.check_reaction(
            "N-alkylation of primary amines with alkyl halides", rsmi
        ) or checker.check_reaction("N-alkylation of secondary amines with alkyl halides", rsmi):
            if "N-alkylation of primary amines with alkyl halides" not in findings_json["atomic_checks"]["named_reactions"]:
                findings_json["atomic_checks"]["named_reactions"].append("N-alkylation of primary amines with alkyl halides")
            if "N-alkylation of secondary amines with alkyl halides" not in findings_json["atomic_checks"]["named_reactions"]:
                findings_json["atomic_checks"]["named_reactions"].append("N-alkylation of secondary amines with alkyl halides")
            print("N-alkylation identified: N-alkylation with alkyl halides")
            return "N-alkylation"

        return "other"

    def dfs_traverse(node, depth=0):
        nonlocal o_alkylation_found, halogen_exchange_found, n_alkylation_found, reaction_sequence

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reaction_type = identify_reaction_type(rsmi)

                print(f"Analyzing reaction at depth {depth}: {rsmi}")
                print(f"Identified as: {reaction_type}")

                # Store reaction type with depth for proper ordering
                reaction_sequence.append((reaction_type, depth))

                if reaction_type == "O-alkylation":
                    o_alkylation_found = True
                elif reaction_type == "halogen-exchange":
                    halogen_exchange_found = True
                elif reaction_type == "N-alkylation":
                    n_alkylation_found = True

        # Process children
        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same (as per new rule)
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth (retrosynthetic order)
    reaction_sequence.sort(key=lambda x: x[1])
    reaction_types = [r[0] for r in reaction_sequence]

    print(f"Reaction sequence: {reaction_types}")
    print(f"O-alkylation found: {o_alkylation_found}")
    print(f"Halogen exchange found: {halogen_exchange_found}")
    print(f"N-alkylation found: {n_alkylation_found}")

    result = False
    # Check if all required reactions are found
    if o_alkylation_found and halogen_exchange_found and n_alkylation_found:
        # Add co-occurrence constraint
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "O-alkylation",
                    "halogen-exchange",
                    "N-alkylation"
                ]
            }
        })

        # Get indices of each reaction type
        try:
            o_alkylation_idx = reaction_types.index("O-alkylation")
            halogen_exchange_idx = reaction_types.index("halogen-exchange")
            n_alkylation_idx = reaction_types.index("N-alkylation")

            # Check the order - in synthetic (forward) direction, we should see:
            # O-alkylation → halogen exchange → N-alkylation
            # But in retrosynthetic direction (as we traverse), we see:
            # N-alkylation → halogen exchange → O-alkylation
            if n_alkylation_idx < halogen_exchange_idx < o_alkylation_idx:
                print("Correct sequence found in retrosynthetic direction")
                result = True
                # Add sequence constraint
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "ordered_targets": [
                            "N-alkylation",
                            "halogen-exchange",
                            "O-alkylation"
                        ],
                        "direction": "retrosynthetic"
                    }
                })
            else:
                print(
                    f"Reactions found but not in the correct sequence: O-alkylation at {o_alkylation_idx}, halogen exchange at {halogen_exchange_idx}, N-alkylation at {n_alkylation_idx}"
                )
        except ValueError:
            # This shouldn't happen since we already checked that all reactions are found
            print("Error: Reaction type missing from sequence despite being found")

    return result, findings_json
