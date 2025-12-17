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
    This function detects the formation of an aldehyde via a Directed Ortho Metalation of Arenes (DoM) strategy, a specific type of C-H functionalization.
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

    ch_to_aldehyde_detected = False
    aldehyde_formed_in_product = False
    reactant_without_aldehyde_present = False
    dom_reaction_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ch_to_aldehyde_detected, aldehyde_formed_in_product, reactant_without_aldehyde_present, dom_reaction_detected, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains aldehyde
                if checker.check_fg("Aldehyde", product_smiles):
                    aldehyde_formed_in_product = True
                    if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")

                    # Check if any reactant does NOT contain aldehyde
                    reactant_without_aldehyde_local = False
                    for reactant in reactants_smiles:
                        if not checker.check_fg("Aldehyde", reactant):
                            reactant_without_aldehyde_local = True
                            break
                    
                    if reactant_without_aldehyde_local:
                        reactant_without_aldehyde_present = True
                        # Check for specific C-H functionalization reactions
                        if checker.check_reaction("Directed ortho metalation of arenes", rsmi):
                            dom_reaction_detected = True
                            if "Directed ortho metalation of arenes" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("Directed ortho metalation of arenes")

                            print(f"Detected C-H functionalization to aldehyde at depth {depth}")
                            print(f"Reaction SMILES: {rsmi}")
                            ch_to_aldehyde_detected = True
                            # No return here, continue to traverse children if needed for other findings

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # After traversal, check for the structural constraint
    if ch_to_aldehyde_detected: # This flag implies all conditions for the strategy were met
        # The strategy JSON defines a co-occurrence constraint for "Directed ortho metalation of arenes" and "aldehyde_formation"
        # The `ch_to_aldehyde_detected` flag already encapsulates this logic.
        # We add the constraint object if the overall strategy is detected.
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Directed ortho metalation of arenes",
                    "aldehyde_formation"
                ],
                "scope": "reaction"
            }
        })

    return ch_to_aldehyde_detected, findings_json
