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


AZIDE_FORMATION_REACTIONS = [
    "Formation of Azides from halogens",
    "Formation of Azides from boronic acids",
    "Amine to azide",
]

AZIDE_UTILIZATION_REACTIONS = [
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Azide to amine reduction (Staudinger)",
    "Azide-nitrile click cycloaddition to tetrazole",
    "Azide-nitrile click cycloaddition to triazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects the use of azide chemistry by identifying specific, well-defined azide formation and utilization reactions. It checks if any reaction in the route matches a predefined list of formation reactions (e.g., from halides, boronic acids, amines) or utilization reactions (e.g., Huisgen cycloadditions, Staudinger reduction, azide-nitrile cycloadditions). A route is flagged if it contains at least one azide intermediate and at least one of these specified formation or utilization reactions.
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

    # Track azide-related transformations
    azide_forming_reactions = 0
    azide_utilizing_reactions = 0
    azide_intermediates = 0

    def dfs_traverse(node, depth=0):
        nonlocal azide_forming_reactions, azide_utilizing_reactions, azide_intermediates, findings_json

        if node["type"] == "mol":
            # Check if molecule contains azide group
            if checker.check_fg("Azide", node["smiles"]):
                azide_intermediates += 1
                if "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Azide")
                print(f"Azide intermediate detected: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for azide formation reactions
            for r_name in AZIDE_FORMATION_REACTIONS:
                if checker.check_reaction(r_name, rsmi):
                    azide_forming_reactions += 1
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    print(f"Azide formation reaction detected: {rsmi}")
                    break # Only count one match per reaction node

            # Check for azide utilization reactions
            for r_name in AZIDE_UTILIZATION_REACTIONS:
                if checker.check_reaction(r_name, rsmi):
                    azide_utilizing_reactions += 1
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    print(f"Azide utilizing reaction detected: {rsmi}")
                    break # Only count one match per reaction node

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., chemical), increase depth
                new_depth = depth + 1
            # If current node is a reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Determine if our strategy is present
    # Strategy is present if we have azide intermediates AND either formation or utilization reactions
    azide_strategy_present = azide_intermediates > 0 and (
        azide_forming_reactions > 0 or azide_utilizing_reactions > 0
    )

    if azide_strategy_present:
        # Add the structural constraint if the overall strategy is detected
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Azide",
                    "azide_formation_or_utilization_reaction"
                ]
            }
        })

    print(f"Azide-based functional group manipulation strategy detected: {azide_strategy_present}")
    print(f"Azide forming reactions: {azide_forming_reactions}")
    print(f"Azide utilizing reactions: {azide_utilizing_reactions}")
    print(f"Azide intermediates: {azide_intermediates}")

    return azide_strategy_present, findings_json
