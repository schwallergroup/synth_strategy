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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy involving a late-stage reductive amination
    to introduce a nitrogen-containing moiety.
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

    has_reductive_amination = False

    def dfs_traverse(node, depth=0):
        nonlocal has_reductive_amination, findings_json

        # Check if this is a reaction node
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Check if this is a late-stage reaction (depth 0 or 1)
            if depth <= 1:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is a reductive amination reaction
                is_reductive_amination_aldehyde = checker.check_reaction(
                    "Reductive amination with aldehyde", rsmi
                )
                is_reductive_amination_ketone = checker.check_reaction(
                    "Reductive amination with ketone", rsmi
                )

                if is_reductive_amination_aldehyde or is_reductive_amination_ketone:
                    print(f"Detected reductive amination at depth {depth}")
                    has_reductive_amination = True
                    if is_reductive_amination_aldehyde:
                        findings_json["atomic_checks"]["named_reactions"].append("Reductive amination with aldehyde")
                    if is_reductive_amination_ketone:
                        findings_json["atomic_checks"]["named_reactions"].append("Reductive amination with ketone")
                    
                    # Add the structural constraint if a late-stage reductive amination is found
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "targets": [
                                "Reductive amination with aldehyde",
                                "Reductive amination with ketone"
                            ],
                            "position": "late_stage",
                            "condition": "any",
                            "qualifier": "depth <= 1"
                        }
                    })

        # Traverse children with incremented depth
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_reductive_amination, findings_json
