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


ESTER_HYDROLYSIS_REACTIONS = [
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
]

ACID_TO_AMIDE_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Carboxylic acid to amide conversion",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving an ester-acid-amide transformation sequence by identifying specific reaction types. It first looks for an ester hydrolysis reaction (e.g., saponification) and subsequently an amide formation from a carboxylic acid (e.g., via coupling agents or acyl chloride intermediates).
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

    # Track transformations with their depths
    transformations = []

    def dfs_traverse(node, depth=0):
        nonlocal transformations, findings_json
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for ester hydrolysis using specific reaction checks
                for rxn_type in ESTER_HYDROLYSIS_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Ester hydrolysis reaction detected at depth {depth}")
                        transformations.append(("ester_hydrolysis", depth))
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)

                # Check for amide formation from an acid (or its activated form like acyl chloride)
                for rxn_type in ACID_TO_AMIDE_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Amide formation reaction detected at depth {depth}")
                        transformations.append(("amide_formation", depth))
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (i.e., it's a chemical), increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Sort transformations by depth (ascending order)
    transformations.sort(key=lambda x: x[1])
    print(f"Sorted transformations: {transformations}")

    # Check if we have both transformations in the correct order for retrosynthesis
    # In retrosynthesis: amide -> acid -> ester (lower depth = later stage in synthesis)
    has_amide_to_acid = False
    has_acid_to_ester = False

    # First, look for amide formation
    for trans, depth in transformations:
        if trans == "amide_formation":
            has_amide_to_acid = True
            # Now look for ester hydrolysis at greater depth (earlier in synthesis)
            for trans2, depth2 in transformations:
                if trans2 == "ester_hydrolysis" and depth2 > depth:
                    has_acid_to_ester = True
                    break
            if has_acid_to_ester:
                break

    has_sequence = has_amide_to_acid and has_acid_to_ester
    if has_sequence:
        print("Ester-acid-amide transformation sequence detected")
        # Add structural constraints if the sequence is detected
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "before": "ester_hydrolysis",
                "after": "amide_formation",
                "description": "An ester hydrolysis reaction must occur at an earlier synthetic step (higher depth) than an amide formation reaction."
            }
        })
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "ester_hydrolysis",
                    "amide_formation"
                ],
                "description": "The route must contain at least one ester hydrolysis reaction and at least one amide formation reaction."
            }
        })
    else:
        print(
            f"Sequence not detected. Amide formation: {has_amide_to_acid}, Ester hydrolysis: {has_acid_to_ester}"
        )

    return has_sequence, findings_json
