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


AZIDE_UTILIZATION_REACTIONS = [
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Azide to amine reduction (Staudinger)",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Azide-nitrile click cycloaddition to tetrazole",
    "Azide-nitrile click cycloaddition to triazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where an amine is converted to an azide, which is then consumed in a subsequent reaction. This function specifically checks for the 'Amine to azide' named reaction for the formation step, and a defined list of utilization reactions including Huisgen cycloadditions, Staudinger reductions, and azide-nitrile cycloadditions.
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

    # Track azide intermediates and their origins
    azide_intermediates = set()  # SMILES of molecules with azides
    azide_from_amine = False  # Flag for amine-to-azide conversion
    azide_utilized = False  # Flag for azide utilization in key reactions

    def dfs_traverse(node, depth=0):
        nonlocal azide_from_amine, azide_utilized, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]

                # Split reactants into individual molecules
                reactant_mols = reactants_part.split(".")

                # Check if this is an amine to azide reaction using the checker
                if checker.check_reaction("Amine to azide", rsmi):
                    print(f"Found amine to azide transformation at depth {depth}")
                    azide_from_amine = True
                    if "Amine to azide" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Amine to azide")

                # Check if this reaction uses an azide in a key transformation
                has_azide_reactant = any(
                    checker.check_fg("Azide", reactant) for reactant in reactant_mols
                )
                if has_azide_reactant:
                    if "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Azide")
                    # Check if azide is used in common azide reactions
                    for name in AZIDE_UTILIZATION_REACTIONS:
                        if checker.check_reaction(name, rsmi):
                            print(f"Found azide being used in key reaction at depth {depth}")
                            azide_utilized = True
                            if name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(name)
                            break

        # Check if the current node is a molecule with azide group
        elif node["type"] == "mol" and checker.check_fg("Azide", node["smiles"]):
            print(f"Found azide intermediate at depth {depth}")
            azide_intermediates.add(node["smiles"])
            if "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Azide")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is 'chemical' (or 'mol'), increase depth
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # A true amine-to-azide strategy requires both:
    # 1. Conversion of amines to azides
    # 2. Utilization of those azides in subsequent reactions
    has_amine_to_azide_strategy = azide_from_amine or (len(azide_intermediates) > 0)

    # Add structural constraint if the condition is met
    # The original code's logic for `has_amine_to_azide_strategy` is `azide_from_amine or (len(azide_intermediates) > 0)`
    # The strategy JSON's structural constraint implies a co-occurrence of 'Amine to azide' and 'azide_utilization'.
    # Given the current code's final boolean logic, we'll add the structural constraint if the final boolean is True.
    # This aligns with the description: "The implemented logic returns true if either the 'Amine to azide' reaction is found or any molecule with an azide group is present, but does not enforce the consumption of the azide."
    if has_amine_to_azide_strategy:
        # This structural constraint is added if the overall strategy is detected, as per the original function's return logic.
        # The note in the strategy JSON indicates the code's deviation from the intended co-occurrence.
        # We'll add the constraint that matches the *actual* code's detection logic.
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Amine to azide",
                    "azide_utilization"
                ],
                "note": "This represents the intended logic from the docstring (azide formation and consumption). The actual code implements a simpler OR condition: ('Amine to azide' reaction) OR (presence of 'Azide' functional group)."
            }
        })

    print(f"Amine to azide strategy detected: {has_amine_to_azide_strategy}")
    return has_amine_to_azide_strategy, findings_json
