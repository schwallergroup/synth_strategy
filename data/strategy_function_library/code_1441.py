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
    Detects a multi-step strategy culminating in a late-stage amide coupling. It identifies routes where the required carboxylic acid is either generated in a prior step (e.g., via ester deprotection or oxidation) or activated to an acyl chloride before the final coupling.
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

    # Track if we found the key elements of the strategy
    found_amide_coupling = False
    found_acid_chloride_formation = False
    found_ester_hydrolysis = False
    found_acid_activation = False  # Track any acid activation method

    # Define the structural constraints from the input JSON for easy access
    structural_constraints_metadata = [
        {
            "type": "positional",
            "details": {
                "event_name": "late_stage_amide_coupling",
                "position": "depth <= 2",
                "target_group": {
                    "logic": "or",
                    "targets": [
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Schotten-Baumann_amide",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analog_N",
                        "Carboxylic acid with primary amine to amide",
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids"
                    ]
                }
            }
        },
        {
            "type": "positional",
            "details": {
                "event_name": "acid_activation",
                "position": "1 <= depth <= 3",
                "target_group": {
                    "logic": "and",
                    "targets": [
                        {
                            "type": "functional_group",
                            "name": "Carboxylic acid",
                            "location": "reactant"
                        },
                        {
                            "type": "functional_group",
                            "name": "Acyl halide",
                            "location": "product"
                        }
                    ]
                }
            }
        },
        {
            "type": "positional",
            "details": {
                "event_name": "acid_generation",
                "position": "2 <= depth <= 4",
                "target_group": {
                    "logic": "or",
                    "targets": [
                        "Ester saponification (alkyl deprotection)",
                        "Ester saponification (methyl deprotection)",
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                        "COOH ethyl deprotection",
                        "Deprotection of carboxylic acid",
                        "Oxidation of aldehydes to carboxylic acids",
                        "Oxidation of nitrile to carboxylic acid"
                    ]
                }
            }
        },
        {
            "type": "co-occurrence",
            "details": {
                "logic": "A and (B or C)",
                "event_mapping": {
                    "A": "late_stage_amide_coupling",
                    "B": "acid_activation",
                    "C": "acid_generation"
                }
            }
        }
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_coupling, found_acid_chloride_formation, found_ester_hydrolysis, found_acid_activation, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                reactants = reactants_str.split(".")

                # Check for amide coupling at late stage (depth 0-2)
                if depth <= 2:
                    amide_coupling_reactions = [
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Schotten-Baumann_amide",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Carboxylic acid with primary amine to amide",
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids"
                    ]
                    for r_name in amide_coupling_reactions:
                        if checker.check_reaction(r_name, rsmi):
                            found_amide_coupling = True
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                            print(f"Found late-stage amide coupling at depth {depth}")
                            # Add the corresponding structural constraint if found
                            if {"type": "positional", "details": {"event_name": "late_stage_amide_coupling", "position": "depth <= 2", "target_group": {"logic": "or", "targets": amide_coupling_reactions}}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(next(item for item in structural_constraints_metadata if item.get("details", {}).get("event_name") == "late_stage_amide_coupling"))
                            break

                # Check for acid chloride formation or other activation (depth 1-3)
                if 1 <= depth <= 3:
                    # Fallback to functional group checking
                    has_carboxylic_acid = False
                    for reactant in reactants:
                        if checker.check_fg("Carboxylic acid", reactant):
                            has_carboxylic_acid = True
                            if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                    has_acid_chloride_product = checker.check_fg("Acyl halide", product_str)
                    if has_acid_chloride_product:
                        if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")

                    if has_carboxylic_acid and has_acid_chloride_product:
                        found_acid_chloride_formation = True
                        found_acid_activation = True
                        print(f"Found acid chloride formation (FG check) at depth {depth}")
                        # Add the corresponding structural constraint if found
                        if {"type": "positional", "details": {"event_name": "acid_activation", "position": "1 <= depth <= 3", "target_group": {"logic": "and", "targets": [{"type": "functional_group", "name": "Carboxylic acid", "location": "reactant"}, {"type": "functional_group", "name": "Acyl halide", "location": "product"}]}}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(next(item for item in structural_constraints_metadata if item.get("details", {}).get("event_name") == "acid_activation"))

                # Check for ester hydrolysis (deprotection) or other acid generation (depth 2-4)
                if 2 <= depth <= 4:
                    # Check for specific reaction types
                    acid_generation_reactions = [
                        "Ester saponification (alkyl deprotection)",
                        "Ester saponification (methyl deprotection)",
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                        "COOH ethyl deprotection",
                        "Deprotection of carboxylic acid",
                        "Oxidation of aldehydes to carboxylic acids",
                        "Oxidation of nitrile to carboxylic acid"
                    ]
                    for r_name in acid_generation_reactions:
                        if checker.check_reaction(r_name, rsmi):
                            found_ester_hydrolysis = True
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                            print(f"Found carboxylic acid generation at depth {depth}")
                            # Add the corresponding structural constraint if found
                            if {"type": "positional", "details": {"event_name": "acid_generation", "position": "2 <= depth <= 4", "target_group": {"logic": "or", "targets": acid_generation_reactions}}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(next(item for item in structural_constraints_metadata if item.get("details", {}).get("event_name") == "acid_generation"))
                            break

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemical nodes)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reaction nodes)
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Print overall results for debugging
    print(f"Amide coupling found: {found_amide_coupling}")
    print(f"Acid chloride formation found: {found_acid_chloride_formation}")
    print(f"Ester hydrolysis/acid generation found: {found_ester_hydrolysis}")

    # Determine the final boolean result
    result = found_amide_coupling and (found_acid_chloride_formation or found_ester_hydrolysis)

    # Add the co-occurrence structural constraint if the overall strategy is found
    if result:
        co_occurrence_constraint = next(item for item in structural_constraints_metadata if item.get("type") == "co-occurrence")
        if co_occurrence_constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(co_occurrence_constraint)

    return result, findings_json
