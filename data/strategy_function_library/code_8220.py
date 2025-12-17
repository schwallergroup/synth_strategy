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


ESTER_FORMATION_REACTIONS = [
    "Esterification of Carboxylic Acids",
    "Schotten-Baumann to ester",
    "O-alkylation of carboxylic acids with diazo compounds",
    "Oxidative esterification of primary alcohols",
    "Transesterification",
    "Acetic anhydride and alcohol to ester",
    "Mitsunobu esterification",
    "Intramolecular transesterification/Lactone formation",
]

ESTER_HYDROLYSIS_REACTIONS = [
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "COOH ethyl deprotection",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Identifies an ester protecting group strategy. This is defined as a sequence where an ester is formed in an early step and subsequently hydrolyzed in a later step. The function specifically checks for formation reactions from the `ESTER_FORMATION_REACTIONS` list and hydrolysis reactions from the `ESTER_HYDROLYSIS_REACTIONS` list, ensuring the formation step occurs at a greater depth (earlier in the synthesis) than the hydrolysis step.
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

    ester_formation_reactions = []
    ester_hydrolysis_reactions = []
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rxn_smiles = node["metadata"]["mapped_reaction_smiles"]

            # Check for ester formation reactions
            for name in ESTER_FORMATION_REACTIONS:
                if checker.check_reaction(name, rxn_smiles):
                    ester_formation_reactions.append((depth, rxn_smiles))
                    if name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                    print(f"Ester formation reaction found at depth {depth}: {rxn_smiles[:50]}...")
                    break # Assuming only one match per reaction node is sufficient

            # Check for ester hydrolysis reactions
            for name in ESTER_HYDROLYSIS_REACTIONS:
                if any(checker.check_reaction(name, rxn_smiles) for name in ESTER_HYDROLYSIS_REACTIONS):
                    ester_hydrolysis_reactions.append((depth, rxn_smiles))
                    if name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                    print(f"Ester hydrolysis reaction found at depth {depth}: {rxn_smiles[:50]}...")
                    break # Assuming only one match per reaction node is sufficient

        # Recursively traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have both ester formation and hydrolysis reactions
    if ester_formation_reactions and ester_hydrolysis_reactions:
        # Add co-occurrence constraint if both types of reactions are found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "any_ester_formation_reaction",
                    "any_ester_hydrolysis_reaction"
                ],
                "description": "Requires the presence of at least one reaction from the ESTER_FORMATION_REACTIONS list and at least one reaction from the ESTER_HYDROLYSIS_REACTIONS list."
            }
        })

        # A protecting group strategy requires formation (high depth) to occur before hydrolysis (low depth).
        # Note: depth=1 is the final step, depth=max_depth is the first step.
        # Therefore, we need to find a formation step with depth > hydrolysis step depth.
        valid_sequence = any(
            f_depth > h_depth
            for f_depth, _ in ester_formation_reactions
            for h_depth, _ in ester_hydrolysis_reactions
        )

        if valid_sequence:
            result = True
            # Add sequence constraint if valid sequence is found
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "any_ester_formation_reaction",
                    "after": "any_ester_hydrolysis_reaction",
                    "description": "An ester formation reaction must occur at a greater depth (earlier in the synthesis) than an ester hydrolysis reaction."
                }
            })
            print(
                f"Ester protecting group strategy detected: {len(ester_formation_reactions)} formation, {len(ester_hydrolysis_reactions)} hydrolysis"
            )
        else:
            print(f"Ester reactions not in correct order for a protecting group strategy.")

    return result, findings_json