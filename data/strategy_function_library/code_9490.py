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


ESTER_HYDROLYSIS_REACTIONS = [
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "COOH ethyl deprotection",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy involving late-stage ester hydrolysis, checked at depths 1 and 2. The strategy is identified if the reaction matches a specific named hydrolysis reaction (defined in `ESTER_HYDROLYSIS_REACTIONS`) or if it results in the net conversion of an ester to a carboxylic acid.
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

    has_late_stage_hydrolysis = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_hydrolysis, findings_json

        if node["type"] == "reaction" and 1 <= depth <= 2:  # Check at depths 1-2 (late stage)
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                print(f"No reaction SMILES found in metadata at depth {depth}")
                return

            print(f"Examining reaction at depth {depth}: {rsmi}")
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Check if this is an ester hydrolysis reaction using known reaction types
                is_hydrolysis = False
                for name in ESTER_HYDROLYSIS_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        is_hydrolysis = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)

                print(f"Is recognized hydrolysis reaction: {is_hydrolysis}")

                # Check for ester in reactants and carboxylic acid in product
                has_ester_reactant = False
                for r in reactants_smiles:
                    if checker.check_fg("Ester", r):
                        has_ester_reactant = True
                        if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ester")
                        break

                has_acid_product = checker.check_fg("Carboxylic acid", product_smiles)
                if has_acid_product:
                    if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                print(f"Has ester in reactants: {has_ester_reactant}")
                print(f"Has carboxylic acid in product: {has_acid_product}")

                # Check if any ester is present in reactants but not in product
                ester_in_product = checker.check_fg("Ester", product_smiles)
                print(f"Has ester in product: {ester_in_product}")

                # Detect hydrolysis either by reaction type or by functional group transformation
                condition_1 = (is_hydrolysis and has_ester_reactant and has_acid_product)
                condition_2 = (has_ester_reactant and has_acid_product and not ester_in_product)

                if condition_1 or condition_2:
                    print(f"Detected late-stage ester hydrolysis at depth {depth}")
                    has_late_stage_hydrolysis = True

                    # Add structural constraints based on which condition was met
                    if 1 <= depth <= 2:
                        positional_constraint = {
                            "type": "positional",
                            "details": {
                                "target": "ester_hydrolysis",
                                "position": "depth_1_or_2"
                            }
                        }
                        if positional_constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(positional_constraint)

                    if condition_1:
                        co_occurrence_constraint_1 = {
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "known_ester_hydrolysis_reaction",
                                    "Ester",
                                    "Carboxylic acid"
                                ],
                                "comment": "This represents the first path to detection: a named reaction with the expected functional groups."
                            }
                        }
                        if co_occurrence_constraint_1 not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(co_occurrence_constraint_1)

                    if condition_2:
                        co_occurrence_constraint_2 = {
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "Ester",
                                    "Carboxylic acid"
                                ],
                                "comment": "This is part of the second path to detection: a net conversion of functional groups."
                            }
                        }
                        if co_occurrence_constraint_2 not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(co_occurrence_constraint_2)

                        negation_constraint = {
                            "type": "negation",
                            "details": {
                                "target": "Ester",
                                "scope": "product",
                                "comment": "This is the other part of the second path, ensuring the ester is consumed."
                            }
                        }
                        if negation_constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(negation_constraint)

            except Exception as e:
                print(f"Error in hydrolysis detection at depth {depth}: {e}")
                traceback.print_exc() # Added for better debugging

        # Process children with incremented depth
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    print("Starting traversal to find late-stage hydrolysis")
    dfs_traverse(route)
    print(f"Final result: {has_late_stage_hydrolysis}")

    return has_late_stage_hydrolysis, findings_json
