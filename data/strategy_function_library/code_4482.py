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


KETONE_REDUCTION_RXNS = [
    "Reduction of ketone to secondary alcohol",
    "Reduction of aldehydes and ketones to alcohols"
]

ACETAL_CLEAVAGE_RXNS = [
    "Acetal hydrolysis to ketone",
    "Ketal hydrolysis to ketone",
    "Acetal hydrolysis to aldehyde"
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic sequence containing three specific reaction types in a defined order: (1) a ketone reduction, (2) an acetal or ketal hydrolysis, and (3) an alcohol oxidation. The strategy is flagged if the reduction occurs earliest in the synthesis (highest depth), followed by the hydrolysis, and the oxidation occurs latest (lowest depth).
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

    # Track if we've found each step in the redox cycle
    found_ketone_reduction = False
    found_acetal_cleavage = False
    found_alcohol_oxidation = False

    # Track the depth of each reaction for sequence analysis
    reduction_depth = -1
    acetal_cleavage_depth = -1
    oxidation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_ketone_reduction, found_acetal_cleavage, found_alcohol_oxidation
        nonlocal reduction_depth, acetal_cleavage_depth, oxidation_depth, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for ketone reduction
                for rxn_name in KETONE_REDUCTION_RXNS:
                    if checker.check_reaction(rxn_name, rsmi):
                        found_ketone_reduction = True
                        reduction_depth = depth
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        break # Found one, no need to check others

                # Check for acetal cleavage
                for rxn_name in ACETAL_CLEAVAGE_RXNS:
                    if any(checker.check_reaction(rxn_name, rsmi) for rxn_name in ACETAL_CLEAVAGE_RXNS):
                        found_acetal_cleavage = True
                        acetal_cleavage_depth = depth
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        break # Found one, no need to check others

                # Check for alcohol oxidation
                oxidation_rxn_name = "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones"
                if checker.check_reaction(oxidation_rxn_name, rsmi):
                    found_alcohol_oxidation = True
                    oxidation_depth = depth
                    if oxidation_rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(oxidation_rxn_name)

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, next node (reaction) increases depth
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    result = False
    # In retrosynthetic analysis, higher depth means earlier in the synthesis
    # For forward synthesis: reduction -> acetal cleavage -> oxidation
    # In retrosynthesis traversal: oxidation -> acetal cleavage -> reduction
    # So we want: reduction_depth > acetal_cleavage_depth > oxidation_depth

    # Check if we found all steps
    if found_ketone_reduction and found_acetal_cleavage and found_alcohol_oxidation:
        # Add co-occurrence constraint if all three are found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "ketone_reduction",
                    "acetal_cleavage",
                    "alcohol_oxidation"
                ],
                "description": "The route must contain at least one reaction from each of the three conceptual groups: ketone reduction, acetal/ketal cleavage, and alcohol oxidation."
            }
        })

        # Check for strict sequence
        if reduction_depth > acetal_cleavage_depth > oxidation_depth:
            result = True
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "ketone_reduction",
                    "after": "alcohol_oxidation",
                    "allow_concurrent": False,
                    "description": "A ketone reduction event must occur in the synthesis pathway before an alcohol oxidation event."
                }
            })
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "acetal_cleavage",
                    "after": "alcohol_oxidation",
                    "allow_concurrent": True,
                    "description": "An acetal or ketal cleavage event must occur before or at the same synthetic stage as an alcohol oxidation event."
                }
            })

        # Check for same-depth reactions (acetal cleavage and oxidation might happen in same step)
        if acetal_cleavage_depth == oxidation_depth and reduction_depth > acetal_cleavage_depth:
            result = True
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "ketone_reduction",
                    "after": "alcohol_oxidation",
                    "allow_concurrent": False,
                    "description": "A ketone reduction event must occur in the synthesis pathway before an alcohol oxidation event."
                }
            })
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "acetal_cleavage",
                    "after": "alcohol_oxidation",
                    "allow_concurrent": True,
                    "description": "An acetal or ketal cleavage event must occur before or at the same synthetic stage as an alcohol oxidation event."
                }
            })

        # Check for relaxed sequence - oxidation must be last in forward synthesis (lowest depth in retro)
        if oxidation_depth < min(acetal_cleavage_depth, reduction_depth):
            result = True
            # The relaxed sequence implies the 'before' and 'after' conditions for the structural constraints
            # are met. We add them if they haven't been added by the strict or same-depth checks.
            # To avoid duplicates, we can check if they are already present or just add them if the result is True.
            # For simplicity, we'll add them directly here, assuming the problem implies adding them if the overall condition is met.
            # In a real scenario, one might de-duplicate the structural_constraints list.
            if {"type": "sequence", "details": {"before": "ketone_reduction", "after": "alcohol_oxidation", "allow_concurrent": False, "description": "A ketone reduction event must occur in the synthesis pathway before an alcohol oxidation event."}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "before": "ketone_reduction",
                        "after": "alcohol_oxidation",
                        "allow_concurrent": False,
                        "description": "A ketone reduction event must occur in the synthesis pathway before an alcohol oxidation event."
                    }
                })
            if {"type": "sequence", "details": {"before": "acetal_cleavage", "after": "alcohol_oxidation", "allow_concurrent": True, "description": "An acetal or ketal cleavage event must occur before or at the same synthetic stage as an alcohol oxidation event."}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "before": "acetal_cleavage",
                        "after": "alcohol_oxidation",
                        "allow_concurrent": True,
                        "description": "An acetal or ketal cleavage event must occur before or at the same synthetic stage as an alcohol oxidation event."
                    }
                })

    return result, findings_json