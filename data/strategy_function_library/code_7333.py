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
    Detects a strategy involving late-stage deprotection of a protected carboxylic acid.
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

    # Track if a late-stage deprotection is found
    late_stage_deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_deprotection_found, findings_json

        # Check if this is a reaction node
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Only consider reactions at low depth (late-stage)
            if depth <= 2:
                # Add positional constraint if met
                if {"type": "positional", "details": {"target": "any_deprotection_reaction", "position": "late_stage", "condition": "depth <= 2"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "any_deprotection_reaction", "position": "late_stage", "condition": "depth <= 2"}})

                # Check if this is a carboxylic acid deprotection reaction
                is_deprotection = False
                deprotection_reactions = [
                    "COOH ethyl deprotection",
                    "Ester saponification (methyl deprotection)",
                    "Ester saponification (alkyl deprotection)",
                    "Deprotection of carboxylic acid"
                ]
                for r_name in deprotection_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        is_deprotection = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                # Extract product
                product = rsmi.split(">")[-1]

                # Check if product has a carboxylic acid group
                has_carboxylic_acid = checker.check_fg("Carboxylic acid", product)
                if has_carboxylic_acid:
                    if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                if is_deprotection and has_carboxylic_acid:
                    late_stage_deprotection_found = True
                    # Add co-occurrence constraint if met
                    if {"type": "co-occurrence", "details": {"targets": ["any_deprotection_reaction", "Carboxylic acid"], "scope": "reaction_produces_functional_group"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["any_deprotection_reaction", "Carboxylic acid"], "scope": "reaction_produces_functional_group"}})

        # Traverse children (in retrosynthetic direction)
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return late_stage_deprotection_found, findings_json
