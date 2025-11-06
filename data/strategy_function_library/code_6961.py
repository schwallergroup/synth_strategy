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
    This function detects a reduction-oxidation sequence where an ester is reduced to an alcohol
    and then the alcohol is oxidized to an aldehyde or ketone.
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

    # Track if we found the pattern
    found_reduction = False
    found_oxidation = False
    reduction_depth = -1
    oxidation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_reduction, found_oxidation, reduction_depth, oxidation_depth, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            # Check for ester reduction to alcohol
            reduction_reaction_name = "Reduction of ester to primary alcohol"
            if checker.check_reaction(reduction_reaction_name, rsmi):
                found_reduction = True
                reduction_depth = depth
                if reduction_reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(reduction_reaction_name)

            # Check for alcohol oxidation to aldehyde/ketone
            oxidation_reaction_name = "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones"
            if checker.check_reaction(oxidation_reaction_name, rsmi):
                found_oxidation = True
                oxidation_depth = depth
                if oxidation_reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(oxidation_reaction_name)

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # In retrosynthesis, oxidation should be found at a lower depth than reduction
    # (meaning it happens after reduction in forward synthesis)
    correct_sequence = found_reduction and found_oxidation and oxidation_depth < reduction_depth

    if found_reduction and found_oxidation:
        # Co-occurrence constraint
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Reduction of ester to primary alcohol",
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones"
                ]
            }
        })

    if correct_sequence:
        # Sequence constraint
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "before": "Reduction of ester to primary alcohol",
                "after": "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones"
            }
        })

    return correct_sequence, findings_json
