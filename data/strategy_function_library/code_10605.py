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
    """Detects a late-stage ester hydrolysis strategy. This is identified by checking for reactions at a depth of 2 or less where an ester functional group in the reactants is converted to a carboxylic acid in the product, confirmed by a named reaction check for hydrolysis."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    late_stage_esterification_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_esterification_found, findings_json

        if (
            node["type"] == "reaction" and depth <= 2
        ):  # Check late-stage reactions (expanded to depth 2)
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Reactants: {reactants}")
                print(f"Product: {product}")

                has_ester_in_reactants = False
                for r in reactants:
                    if checker.check_fg("Ester", r):
                        has_ester_in_reactants = True
                        if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ester")
                        break

                has_acid_in_product = False
                if checker.check_fg("Carboxylic acid", product):
                    has_acid_in_product = True
                    if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                print(f"Has ester in reactants: {has_ester_in_reactants}")
                print(f"Has carboxylic acid in product: {has_acid_in_product}")

                is_ester_hydrolysis = False
                hydrolysis_reactions = [
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    "Ester saponification (methyl deprotection)",
                    "Ester saponification (alkyl deprotection)",
                    "COOH ethyl deprotection"
                ]
                detected_hydrolysis_reaction = None

                for r_name in hydrolysis_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        is_ester_hydrolysis = True
                        detected_hydrolysis_reaction = r_name
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        break

                print(f"Is ester hydrolysis reaction: {is_ester_hydrolysis}")

                # Detect late-stage ester hydrolysis strategy
                if (
                    has_ester_in_reactants
                    and has_acid_in_product
                    and is_ester_hydrolysis
                ):
                    print(f"Found late-stage ester hydrolysis strategy at depth {depth}")
                    late_stage_esterification_found = True

                    # Add structural constraints
                    if depth <= 2:
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "ester_hydrolysis_reaction",
                                "position": "depth <= 2"
                            }
                        })
                    
                    # Add co-occurrence constraint based on the detected reaction
                    if detected_hydrolysis_reaction:
                        constraint_details = {
                            "scope": "reaction_step",
                            "targets": [
                                "Ester",
                                "Carboxylic acid",
                                detected_hydrolysis_reaction
                            ]
                        }
                        # Check if this specific constraint already exists to avoid duplicates
                        if constraint_details not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": constraint_details
                            })

            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        # Process children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": 
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal
    print("Starting traversal of synthesis route")
    dfs_traverse(route)
    print(f"Late-stage esterification detected: {late_stage_esterification_found}")

    return late_stage_esterification_found, findings_json
