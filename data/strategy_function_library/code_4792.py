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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage reductive amination involving an aromatic amine.
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

    found_strategy = False
    reductive_amination_found = False
    aromatic_amine_found = False
    is_late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal found_strategy, reductive_amination_found, aromatic_amine_found, is_late_stage, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            product_part = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction (depth <= 1)
            if depth <= 1:
                is_late_stage = True
                # Check if this is a reductive amination reaction using the checker
                current_is_reductive_amination = False
                ra_reactions = [
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Reductive amination with alcohol"
                ]
                for ra_name in ra_reactions:
                    if checker.check_reaction(ra_name, rsmi):
                        current_is_reductive_amination = True
                        reductive_amination_found = True
                        if ra_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(ra_name)

                if current_is_reductive_amination:
                    # Check if an aromatic amine (aniline-type) is present in the product.
                    # For a reductive amination, this implies an aromatic amine was a reactant.
                    product_has_aromatic_amine = checker.check_fg("Aniline", product_part)

                    if product_has_aromatic_amine:
                        aromatic_amine_found = True
                        if "Aniline" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aniline")

                        found_strategy = True

        # Continue DFS traversal
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Add structural constraints based on the flags
    if reductive_amination_found and is_late_stage:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": [
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Reductive amination with alcohol"
                ],
                "position": "late_stage"
            }
        })

    if reductive_amination_found and aromatic_amine_found:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    [
                        "Reductive amination with aldehyde",
                        "Reductive amination with ketone",
                        "Reductive amination with alcohol"
                    ],
                    "Aniline"
                ]
            }
        })

    return found_strategy, findings_json
