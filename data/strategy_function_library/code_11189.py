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


AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with ammonia to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage amide formation strategy by checking if the final synthetic step is a recognized amide-forming reaction.
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

    # Track if we found late-stage amide formation
    found_late_stage_amide = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_amide, findings_json

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate reaction
            try:
                # Extract reactants and product
                if "rsmi" not in node.get("metadata", {}):
                    print(f"No rsmi in metadata for reaction node at depth {depth}")
                    return

                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Processing reaction at depth {depth}: {rsmi}")

                # Check if this is an amide formation reaction using the checker function
                is_amide_formation = False
                for rxn_type in AMIDE_FORMATION_REACTIONS:
                    check_result = checker.check_reaction(rxn_type, rsmi)
                    print(f"Checking reaction type {rxn_type}: {check_result}")
                    if check_result:
                        print(f"Matched reaction type: {rxn_type}")
                        is_amide_formation = True
                        # Record the found reaction
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if is_amide_formation:
                    found_late_stage_amide = True
                    print(
                        f"Detected late-stage amide formation in reaction at depth {depth}: {rsmi}"
                    )
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            # Depth remains the same when traversing from reaction to chemical
            new_depth = depth
            if node['type'] == 'chemical':
                new_depth = depth + 1
            # If node['type'] is 'reaction', new_depth remains 'depth'
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    print("Starting traversal of synthesis route")
    dfs_traverse(route)
    print(f"Late-stage amide formation detected: {found_late_stage_amide}")

    # Add structural constraint if late-stage amide formation was found
    if found_late_stage_amide:
        # This corresponds to the 'positional' structural constraint in the strategy JSON
        # We need to ensure this is added only once and matches the definition.
        # The specific details of the constraint are hardcoded here based on the provided JSON.
        structural_constraint_to_add = {
            "type": "positional",
            "details": {
                "targets": [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with ammonia to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with ammonia to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines"
                ],
                "position": "last_or_penultimate_stage"
            }
        }
        if structural_constraint_to_add not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(structural_constraint_to_add)

    return found_late_stage_amide, findings_json
