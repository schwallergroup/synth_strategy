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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage amide bond formation from an ester precursor.
    This strategy is identified in the retrosynthetic direction by finding a
    reaction where an amide in the target molecule is disconnected to an ester,
    and the reaction is classified as an aminolysis or a related ester-to-amide
    conversion.
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

    found_conversion = False

    def dfs_traverse(node, depth=0):
        nonlocal found_conversion, findings_json

        if node["type"] == "reaction" and depth <= 2:  # Focus on late-stage reactions
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if any reactant has an ester group (retrosynthetic perspective)
                reactant_with_ester = None
                for r_smiles in reactants_smiles:
                    if checker.check_fg("Ester", r_smiles):
                        reactant_with_ester = r_smiles
                        if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ester")
                        print(f"Found ester in reactant: {r_smiles}")
                        break

                # Check if product has an amide group (retrosynthetic perspective)
                product_has_amide = False
                if checker.check_fg("Primary amide", product_smiles):
                    product_has_amide = True
                    if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                if checker.check_fg("Secondary amide", product_smiles):
                    product_has_amide = True
                    if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                if checker.check_fg("Tertiary amide", product_smiles):
                    product_has_amide = True
                    if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                if product_has_amide:
                    print(f"Found amide in product: {product_smiles}")

                # Check for relevant reaction types (in retrosynthetic direction)
                is_aminolysis = checker.check_reaction("Aminolysis of esters", rsmi)
                is_ester_to_amide = False
                
                if is_aminolysis:
                    if "Aminolysis of esters" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Aminolysis of esters")
                    if {"type": "positional", "details": {"target": "Aminolysis of esters", "position": "late_stage (depth <= 2)"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Aminolysis of esters", "position": "late_stage (depth <= 2)"}})

                if checker.check_reaction("Ester with ammonia to amide", rsmi):
                    is_ester_to_amide = True
                    if "Ester with ammonia to amide" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Ester with ammonia to amide")
                    if {"type": "positional", "details": {"target": "Ester with ammonia to amide", "position": "late_stage (depth <= 2)"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Ester with ammonia to amide", "position": "late_stage (depth <= 2)"}})

                if checker.check_reaction("Ester with primary amine to amide", rsmi):
                    is_ester_to_amide = True
                    if "Ester with primary amine to amide" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Ester with primary amine to amide")
                    if {"type": "positional", "details": {"target": "Ester with primary amine to amide", "position": "late_stage (depth <= 2)"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Ester with primary amine to amide", "position": "late_stage (depth <= 2)"}})

                if checker.check_reaction("Ester with secondary amine to amide", rsmi):
                    is_ester_to_amide = True
                    if "Ester with secondary amine to amide" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Ester with secondary amine to amide")
                    if {"type": "positional", "details": {"target": "Ester with secondary amine to amide", "position": "late_stage (depth <= 2)"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Ester with secondary amine to amide", "position": "late_stage (depth <= 2)"}})

                if (
                    reactant_with_ester
                    and product_has_amide
                    and (is_aminolysis or is_ester_to_amide)
                ):
                    # This is a retrosynthetic step that disconnects an amide to an ester
                    # Which means in forward synthesis, it's an amide to ester conversion strategy
                    print(f"Confirmed ester to amide conversion (retrosynthetic) at depth {depth}")
                    found_conversion = True
                    # Add the co-occurrence structural constraint
                    co_occurrence_constraint = {
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "Ester",
                                "Primary amide",
                                "Secondary amide",
                                "Tertiary amide",
                                "Aminolysis of esters",
                                "Ester with ammonia to amide",
                                "Ester with primary amine to amide",
                                "Ester with secondary amine to amide"
                            ],
                            "context": "A single reaction step must involve an Ester as a reactant, an Amide (Primary, Secondary, or Tertiary) as a product, and be one of the specified reaction types."
                        }
                    }
                    if co_occurrence_constraint not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(co_occurrence_constraint)

            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return found_conversion, findings_json
