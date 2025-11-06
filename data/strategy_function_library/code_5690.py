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


ALCOHOL_TO_CHLORIDE_REACTIONS = [
    "Alcohol to chloride_sulfonyl chloride",
    "Alcohol to chloride_SOCl2",
    "Alcohol to chloride_CHCl3",
    "Alcohol to chloride_CH2Cl2",
    "Alcohol to chloride_PCl5_ortho",
    "Alcohol to chloride_POCl3_ortho",
    "Alcohol to chloride_POCl3_para",
    "Alcohol to chloride_POCl3",
    "Alcohol to chloride_HCl",
    "Alcohol to chloride_Salt",
    "Alcohol to chloride_Other",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage functional group interconversion from a hydroxyl group to a chloride. This is identified by checking if the reaction is classified as one of the specific reaction types defined in the ALCOHOL_TO_CHLORIDE_REACTIONS list.
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

    # Track if we've seen hydroxyl to chloride conversion
    hydroxyl_to_chloride_seen = False

    # Track if it happened at a late stage (depth <= 1)
    late_stage_conversion = False

    def dfs_traverse(node, depth):
        nonlocal hydroxyl_to_chloride_seen, late_stage_conversion, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and product from reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for hydroxyl to chloride conversion using specific reaction types
                is_alcohol_to_chloride = False
                for rxn_type in ALCOHOL_TO_CHLORIDE_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_alcohol_to_chloride = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)

                # Detect hydroxyl to chloride conversion
                if is_alcohol_to_chloride:
                    print(f"Detected hydroxyl to chloride conversion at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    hydroxyl_to_chloride_seen = True

                    # Check if it's late stage (depth <= 1)
                    if depth <= 1:
                        print(f"This is a late-stage conversion (depth={depth})")
                        late_stage_conversion = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Determine the new depth for recursive calls
        new_depth = depth
        if node["type"] != "reaction": # This means it's a chemical node
            new_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route, 0)

    print(f"Hydroxyl to chloride conversion seen: {hydroxyl_to_chloride_seen}")
    print(f"Late stage conversion: {late_stage_conversion}")

    # Add structural constraint if both conditions are met
    if hydroxyl_to_chloride_seen and late_stage_conversion:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Alcohol to chloride conversion",
                "position": "last_two_stages (depth <= 1)"
            }
        })

    # Return True if we've seen a late-stage hydroxyl to chloride conversion
    result = hydroxyl_to_chloride_seen and late_stage_conversion
    return result, findings_json
