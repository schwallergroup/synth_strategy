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
    This function detects a synthetic strategy where a sulfonamide group is formed
    in the late stage of the synthesis (low depth) from an amine precursor.
    """
    sulfonamide_formation_detected = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Strategy JSON for reference to populate structural_constraints
    strategy_json = {
        "function_id": "code_6775",
        "filepath": "../data/merged_good_perf/code_6775.py",
        "description": "This function detects a synthetic strategy where a sulfonamide group is formed in the late stage of the synthesis (low depth) from an amine precursor.",
        "atomic_checks": {
            "named_reactions": [
                "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                "Sulfonamide synthesis (Schotten-Baumann) secondary amine"
            ],
            "ring_systems": [],
            "functional_groups": [
                "Sulfonamide"
            ]
        },
        "structural_constraints": [
            {
                "type": "positional",
                "details": {
                    "target": "sulfonamide_formation",
                    "position": "late_stage",
                    "condition": "depth <= 2"
                }
            },
            {
                "type": "negation",
                "details": {
                    "target": "Sulfonamide",
                    "scope": "reactants",
                    "context": "The reactants of the sulfonamide formation step must not contain a sulfonamide group."
                }
            },
            {
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "Sulfonamide",
                        "Schotten-Baumann sulfonamide synthesis"
                    ],
                    "scope": "single_reaction",
                    "relationship": "The reaction must be a named Schotten-Baumann synthesis and its product must contain a Sulfonamide functional group."
                }
            }
        ]
    }

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_formation_detected, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Use traversal depth if metadata depth is not available
            node_depth = depth
            if (
                "metadata" in node
                and "depth" in node["metadata"]
                and node["metadata"]["depth"] is not None
            ):
                node_depth = node["metadata"]["depth"]

            # Consider reactions up to 2 steps from final product as late stage
            if node_depth <= 2:
                try:
                    # Check if product contains sulfonamide
                    has_sulfonamide_product = checker.check_fg("Sulfonamide", product_smiles)

                    if has_sulfonamide_product:
                        if "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")

                        # Check if any reactant already has sulfonamide
                        reactant_has_sulfonamide = any(
                            checker.check_fg("Sulfonamide", r) for r in reactants_smiles
                        )

                        # Only proceed if sulfonamide is formed in this reaction
                        if not reactant_has_sulfonamide:
                            # Add negation constraint if this condition is met
                            negation_constraint = next((c for c in strategy_json["structural_constraints"] if c["type"] == "negation" and c["details"]["target"] == "Sulfonamide"), None)
                            if negation_constraint and negation_constraint not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(copy.deepcopy(negation_constraint))

                            # Check if this is a sulfonamide formation reaction
                            is_primary_sulfonamide_rxn = checker.check_reaction(
                                "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                            )
                            is_secondary_sulfonamide_rxn = checker.check_reaction(
                                "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                            )

                            if is_primary_sulfonamide_rxn or is_secondary_sulfonamide_rxn:
                                sulfonamide_formation_detected = True

                                # Add named reactions to findings
                                if is_primary_sulfonamide_rxn and "Sulfonamide synthesis (Schotten-Baumann) primary amine" not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append("Sulfonamide synthesis (Schotten-Baumann) primary amine")
                                if is_secondary_sulfonamide_rxn and "Sulfonamide synthesis (Schotten-Baumann) secondary amine" not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append("Sulfonamide synthesis (Schotten-Baumann) secondary amine")

                                # Add positional constraint if depth condition is met
                                positional_constraint = next((c for c in strategy_json["structural_constraints"] if c["type"] == "positional" and c["details"]["target"] == "sulfonamide_formation"), None)
                                if positional_constraint and positional_constraint not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append(copy.deepcopy(positional_constraint))

                                # Add co-occurrence constraint
                                co_occurrence_constraint = next((c for c in strategy_json["structural_constraints"] if c["type"] == "co-occurrence" and "Sulfonamide" in c["details"]["targets"]), None)
                                if co_occurrence_constraint and co_occurrence_constraint not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append(copy.deepcopy(co_occurrence_constraint))

                except Exception as e:
                    print(f"Error processing reaction SMILES: {str(e)}")

        # Process children with incremented depth
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # Depth increases only when going from chemical to reaction
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return sulfonamide_formation_detected, findings_json
