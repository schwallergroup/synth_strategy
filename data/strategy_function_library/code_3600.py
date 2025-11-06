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


ACYLATION_REACTIONS = [
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Schotten-Baumann_amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
]

NITRATION_REACTIONS = [
    "Aromatic nitration with HNO3",
    "Aromatic nitration with NO3 salt",
    "Aromatic nitration with NO2 salt",
    "Aromatic nitration with alkyl NO2",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a protecting group strategy where an amine is acylated before a subsequent aromatic nitration.
    This function verifies that the acylation step occurs earlier in the synthetic route than the nitration step.
    It identifies the relevant reactions by checking against the predefined named reaction lists `ACYLATION_REACTIONS` and `NITRATION_REACTIONS`.
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

    has_amine_acylation = False
    has_aromatic_nitration = False
    correct_sequence = False
    acylation_depth = -1
    nitration_depth = -1

    def dfs_traverse(node, current_depth=0):
        nonlocal has_amine_acylation, has_aromatic_nitration, correct_sequence
        nonlocal acylation_depth, nitration_depth, findings_json

        if node.get("type") == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            print(f"Analyzing reaction at depth {current_depth}: {rsmi}")

            # Check for acylation reaction
            is_acylation = False
            for reaction_type in ACYLATION_REACTIONS:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Detected acylation reaction: {reaction_type}")
                    is_acylation = True
                    if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    break

            if is_acylation:
                has_amine_acylation = True
                acylation_depth = current_depth
                print(f"Confirmed amine acylation at depth {current_depth}")

            # Check for aromatic nitration using reaction checker
            is_nitration = False
            for reaction_type in NITRATION_REACTIONS:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Detected nitration reaction: {reaction_type}")
                    is_nitration = True
                    if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    break

            if is_nitration:
                has_aromatic_nitration = True
                nitration_depth = current_depth
                print(f"Confirmed aromatic nitration at depth {current_depth}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            # New logic for depth calculation
            if node.get("type") == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, current_depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have both transformations in the correct sequence
    if has_amine_acylation and has_aromatic_nitration:
        if nitration_depth < acylation_depth:
            correct_sequence = True
            print(
                f"Confirmed correct sequence: acylation at depth {acylation_depth}, nitration at depth {nitration_depth}"
            )
            # Add the structural constraint if the sequence is correct
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "amine_acylation",
                    "after": "aromatic_nitration"
                }
            })
        else:
            print(
                f"Incorrect sequence: acylation at depth {acylation_depth}, nitration at depth {nitration_depth}"
            )
    else:
        if not has_amine_acylation:
            print("Missing amine acylation step")
        if not has_aromatic_nitration:
            print("Missing aromatic nitration step")

    return correct_sequence, findings_json
