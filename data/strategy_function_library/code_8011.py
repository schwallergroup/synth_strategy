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


import re

# Refactoring for Enumeration: Isolate the list of reactions.
AROMATIC_FUNCTIONALIZATIONS = [
    "Aromatic nitration with HNO3",
    "Aromatic nitration with NO3 salt",
    "Aromatic nitration with NO2 salt",
    "Aromatic nitration with alkyl NO2",
    "Aromatic chlorination",
    "Aromatic bromination",
    "Aromatic iodination",
    "Aromatic fluorination",
    "Friedel-Crafts acylation",
    "Friedel-Crafts alkylation",
    "Reduction of nitro groups to amines",
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
    "Aromatic hydroxylation",
    "Aromatic sulfonyl chlorination",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis involves at least two sequential aromatic functionalization steps. The identification is based on a predefined list of named reactions, including nitration, halogenation, Friedel-Crafts acylation/alkylation, and others.
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

    # Track functionalization steps
    functionalization_steps = []

    def dfs_traverse(node, depth=0):
        nonlocal functionalization_steps, findings_json
        if node["type"] == "reaction":
            try:
                # Extract reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for aromatic functionalization reactions
                for reaction_type in AROMATIC_FUNCTIONALIZATIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        functionalization_steps.append((reaction_type, depth))
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        print(f"Detected {reaction_type} at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is chemical, depth increases for reaction children
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Sort by depth (higher depth = earlier in synthesis)
    functionalization_steps.sort(key=lambda x: x[1], reverse=True)

    # Check if we have at least 2 sequential aromatic functionalizations
    has_strategy = len(functionalization_steps) >= 2

    if has_strategy:
        step_types = [step[0] for step in functionalization_steps]
        print(f"Aromatic functionalization sequence: {' → '.join(step_types)}")
        # Add the structural constraint if the strategy is found
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "aromatic_functionalization_reaction",
                "operator": ">=",
                "value": 2
            }
        })

    return has_strategy, findings_json
