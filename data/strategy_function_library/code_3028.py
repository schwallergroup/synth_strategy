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
    This function detects a synthetic strategy involving a sequence of nitrogen
    functional group interconversions (nitro → amine → secondary amine).
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

    # Track functional group transformations
    transformations = []
    result = False

    def dfs_traverse(node, path=None):
        nonlocal transformations, findings_json
        if path is None:
            path = []

        current_path = path.copy()

        if node["type"] == "mol":
            current_path.append(node["smiles"])

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]
            depth = node.get("metadata", {}).get("depth", -1)

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # In retrosynthesis, product is the starting material and reactants are the results
            # Check for nitro reduction (nitro → amine)
            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                print(f"Detected nitro reduction reaction at depth {depth}")
                transformations.append(("nitro_to_amine", depth))
                findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

            # Check for primary amine alkylation (primary amine → secondary amine)
            # Include reductive amination which also converts primary to secondary amines
            elif checker.check_reaction(
                "N-alkylation of primary amines with alkyl halides", rsmi
            ):
                print(f"Detected primary amine alkylation/reductive amination at depth {depth}")
                transformations.append(("amine_to_secondary_amine", depth))
                findings_json["atomic_checks"]["named_reactions"].append("N-alkylation of primary amines with alkyl halides")
            elif checker.check_reaction("Reductive amination with aldehyde", rsmi):
                print(f"Detected primary amine alkylation/reductive amination at depth {depth}")
                transformations.append(("amine_to_secondary_amine", depth))
                findings_json["atomic_checks"]["named_reactions"].append("Reductive amination with aldehyde")

        # Traverse children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                new_depth = node.get("metadata", {}).get("depth", -1)
            else:
                # Depth increases when traversing from chemical to reaction
                new_depth = node.get("metadata", {}).get("depth", -1) + 1
            
            # Update the child's metadata with the new depth before recursive call
            if "metadata" not in child:
                child["metadata"] = {}
            child["metadata"]["depth"] = new_depth

            dfs_traverse(child, current_path)

    # Start traversal from root
    # Initialize root depth to 0 or -1 as per original logic if not present
    if "metadata" not in route:
        route["metadata"] = {}
    if "depth" not in route["metadata"]:
        route["metadata"]["depth"] = 0 # Or -1, depending on desired initial depth

    dfs_traverse(route)

    # Sort transformations by depth (highest to lowest, as higher depth is earlier in synthesis)
    # If depth is -1, treat it as a very high value to place it at the beginning (earliest)
    transformations.sort(key=lambda x: -float("inf") if x[1] == -1 else -x[1])

    print(f"All transformations (sorted by depth): {transformations}")

    # Check if we have the nitro → amine → secondary amine sequence
    # In retrosynthesis, this would appear as secondary_amine → amine → nitro
    transformation_types = [t[0] for t in transformations]

    if (
        "nitro_to_amine" in transformation_types
        and "amine_to_secondary_amine" in transformation_types
    ):
        nitro_idx = transformation_types.index("nitro_to_amine")
        amine_idx = transformation_types.index("amine_to_secondary_amine")

        # Check if the transformations are in the correct order
        # In retrosynthesis, amine_to_secondary_amine should come before nitro_to_amine
        # This means amine_idx should be less than nitro_idx in the sorted list
        if amine_idx < nitro_idx:
            print(
                "Detected nitrogen functional group interconversion sequence: nitro → amine → secondary amine"
            )
            result = True
            # Add the structural constraint to findings_json
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "description": "Checks for a sequence where a primary amine is converted to a secondary amine before a nitro group is reduced to that primary amine in the retrosynthetic path.",
                    "before": {
                        "any_of": [
                            "N-alkylation of primary amines with alkyl halides",
                            "Reductive amination with aldehyde"
                        ]
                    },
                    "after": {
                        "all_of": [
                            "Reduction of nitro groups to amines"
                        ]
                    }
                }
            })

    return result, findings_json
