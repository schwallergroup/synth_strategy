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


REDUCTION_REACTIONS_FOR_PIPERAZINE_SEQ = [
    "Reduction of aldehydes and ketones to alcohols",
    "Reduction of carboxylic acid to primary alcohol",
    "Reduction of ester to primary alcohol",
    "Reduction of ketone to secondary alcohol",
    "Reduction of nitrile to amine",
    "Hydrogenation (double to single)",
    "Hydrogenation (triple to double)",
    "Arene hydrogenation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthesis route involving a piperazine scaffold contains a reduction reaction followed by a dehydration reaction. The reduction reactions checked are: Reduction of aldehydes and ketones to alcohols, Reduction of carboxylic acid to primary alcohol, Reduction of ester to primary alcohol, Reduction of ketone to secondary alcohol, Reduction of nitrile to amine, Hydrogenation (double to single), Hydrogenation (triple to double), and Arene hydrogenation. Dehydration is identified as the conversion of an alcohol to an alkene.
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

    has_piperazine = False
    has_reduction = False
    has_dehydration = False
    reduction_depth = -1
    dehydration_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal has_piperazine, has_reduction, has_dehydration, reduction_depth, dehydration_depth, findings_json

        if node["type"] == "mol":
            if not has_piperazine and checker.check_ring("piperazine", node["smiles"]):
                has_piperazine = True
                if "piperazine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("piperazine")

        elif node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            # Check for reduction reactions
            for rxn in REDUCTION_REACTIONS_FOR_PIPERAZINE_SEQ:
                if checker.check_reaction(rxn, rsmi):
                    has_reduction = True
                    reduction_depth = depth
                    if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                    break # Found one reduction, no need to check others for this reaction node

            # Check for dehydration (alcohol elimination to alkene)
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]
            
            alcohol_in_reactants = False
            for r in reactants:
                if checker.check_fg("Primary alcohol", r):
                    alcohol_in_reactants = True
                    if "Primary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                if checker.check_fg("Secondary alcohol", r):
                    alcohol_in_reactants = True
                    if "Secondary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
                if checker.check_fg("Tertiary alcohol", r):
                    alcohol_in_reactants = True
                    if "Tertiary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")

            alkene_in_product = checker.check_fg("Alkene", product)
            if alkene_in_product:
                if "Alkene" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Alkene")

            if alcohol_in_reactants and alkene_in_product:
                has_dehydration = True
                dehydration_depth = depth
                if "dehydration" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("dehydration")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # In forward synthesis, reduction must happen before dehydration.
    # A smaller depth value means a later step in the synthesis.
    # So, dehydration_depth should be smaller than reduction_depth.
    correct_sequence = has_reduction and has_dehydration and dehydration_depth < reduction_depth

    result = has_piperazine and correct_sequence

    if has_piperazine and has_reduction and has_dehydration:
        # Co-occurrence constraint
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "piperazine",
                    "any_reduction_from_list",
                    "dehydration"
                ]
            }
        })

    if has_reduction and has_dehydration and dehydration_depth < reduction_depth:
        # Sequence constraint
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "before": "any_reduction_from_list",
                "after": "dehydration"
            }
        })

    return result, findings_json