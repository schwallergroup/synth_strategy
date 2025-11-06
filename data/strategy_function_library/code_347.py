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


CROSS_COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters OTf",
    "Negishi coupling",
    "Stille reaction_aryl",
    "Stille reaction_vinyl",
    "Stille reaction_benzyl",
    "Stille reaction_allyl",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Aryllithium cross-coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthesis involves a late-stage cross-coupling reaction from a predefined list of named reactions. The checked reactions include various types of Suzuki, Negishi, Stille, Hiyama-Denmark, Kumada, and Aryllithium cross-couplings. A reaction is considered 'late-stage' if it occurs within the final two steps of the synthesis (depth 0 or 1).
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

    has_late_stage_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_coupling, findings_json

        if node["type"] == "reaction" and depth <= 1:
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for a specific list of named cross-coupling reactions.
                if not has_late_stage_coupling:
                    # print(f"Checking reaction: {rsmi}") # Original print statement removed
                    for rxn in CROSS_COUPLING_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            has_late_stage_coupling = True
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            # Add the structural constraint here as it's tied to the detection of a late-stage coupling
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "targets": [
                                        "Suzuki coupling with boronic acids",
                                        "Suzuki coupling with boronic esters",
                                        "Suzuki coupling with boronic acids OTf",
                                        "Suzuki coupling with boronic esters OTf",
                                        "Negishi coupling",
                                        "Stille reaction_aryl",
                                        "Stille reaction_vinyl",
                                        "Stille reaction_benzyl",
                                        "Stille reaction_allyl",
                                        "Hiyama-Denmark Coupling",
                                        "Kumada cross-coupling",
                                        "Aryllithium cross-coupling"
                                    ],
                                    "position": "final_two_stages"
                                }
                            })
                            break # Exit loop once a match is found and recorded

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth if node["type"] == "mol" else depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return has_late_stage_coupling, findings_json
