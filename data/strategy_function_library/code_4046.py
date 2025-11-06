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
    This function detects if the synthetic route incorporates a piperazine moiety
    in a late stage of the synthesis.
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

    # Track if we found a late-stage piperazine incorporation
    found_piperazine_incorporation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_piperazine_incorporation, findings_json

        # A reaction node at depth 1, 2, or 3 is considered "late-stage".
        # Depth=1 is the final reaction step.
        if node.get("type") == "reaction" and 1 <= depth <= 3:
            # Assumes the node object contains structured product/reactant info.
            product = node.get("products", [None])[0]
            reactants = node.get("reactants", [])

            if not product or not reactants:
                return # Cannot analyze incomplete reaction data

            # Check if product contains piperazine
            if checker.check_ring("piperazine", product):
                findings_json["atomic_checks"]["ring_systems"].append("piperazine")
                # "Incorporation" means the product has a piperazine ring that was not
                # present in ALL reactants. This correctly identifies de novo formation
                # and coupling of a piperazine-containing fragment, while excluding
                # simple modifications to a molecule that already contains piperazine.
                if not all(checker.check_ring("piperazine", r) for r in reactants):
                    found_piperazine_incorporation = True
                    # This implies a ring formation or incorporation, which aligns with 'ring_formation'
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                    
                    # Add the structural constraint if the condition is met
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "ring_formation",
                            "position": "late_stage"
                        }
                    })

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node.get("type") != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return found_piperazine_incorporation, findings_json
