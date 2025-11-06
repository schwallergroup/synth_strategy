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
    This function detects if the synthesis includes a late-stage amidation, defined as the formation of an amide from a carboxylic acid or its derivative (e.g., acyl chloride, ester) in the last or second-to-last step.
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

    has_late_amidation = False

    def dfs_traverse(node, depth=0, path=None):
        nonlocal has_late_amidation, findings_json

        if path is None:
            path = []

        # Add current node to path
        path.append(node)

        # Check if this is a reaction node at depth 1 or 2 (late-stage)
        if node["type"] == "reaction" and 1 <= depth <= 2:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is an amidation reaction using comprehensive list
                amidation_reactions = [
                    "Carboxylic acid with primary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Ester with primary amine to amide",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with secondary amine to amide",
                    "Acyl chloride with ammonia to amide",
                    "Ester with ammonia to amide",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                ]

                # Define functional groups to check
                acid_groups = ["Carboxylic acid", "Acyl halide", "Ester", "Anhydride"]
                amide_groups = ["Primary amide", "Secondary amide", "Tertiary amide"]

                is_amidation = False
                has_acid = False
                has_amide = False

                # Check if any of the known amidation reactions match
                for reaction in amidation_reactions:
                    if checker.check_reaction(reaction, rsmi):
                        print(f"Matched reaction type: {reaction}")
                        is_amidation = True
                        if reaction not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction)

                # Check for acid groups in reactants
                for reactant in reactants:
                    for acid in acid_groups:
                        if checker.check_fg(acid, reactant):
                            print(f"Found {acid} in reactant: {reactant}")
                            has_acid = True
                            if acid not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(acid)

                # Check for amide groups in product
                for amide in amide_groups:
                    if checker.check_fg(amide, product):
                        print(f"Found {amide} in product: {product}")
                        has_amide = True
                        if amide not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(amide)

                # If we didn't match a known reaction but found acid → amide transformation
                if not is_amidation and has_acid and has_amide:
                    print("Detected amidation based on functional group transformation")
                    is_amidation = True

                if is_amidation and has_acid and has_amide:
                    has_late_amidation = True
                    print(f"Confirmed late-stage amidation at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    # Add the structural constraint if the condition is met
                    if {"type": "positional", "details": {"target": "amidation", "position": "last_or_second_to_last_stage"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amidation", "position": "last_or_second_to_last_stage"}})

        # Continue traversing
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same

            # Create a new path for each branch to avoid modifying the current path
            dfs_traverse(child, new_depth, path.copy())

    dfs_traverse(route)
    print(
        f"Synthesis {'includes' if has_late_amidation else 'does not include'} late-stage amidation"
    )
    return has_late_amidation, findings_json
