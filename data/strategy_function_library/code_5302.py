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


AMINE_FUNCTIONALIZATION_REACTIONS = [
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving nitro reduction followed by
    amine functionalization (e.g., sulfonamide formation).
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

    # Track reactions and their depths
    reactions_by_depth = {}
    result = False

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Store reaction info by depth
                reactions_by_depth[depth] = {
                    "reactants": reactants,
                    "product": product,
                    "rsmi": rsmi,
                }
            except Exception as e:
                print(f"Error extracting reaction data: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a chemical node
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Sort depths to analyze reaction sequence
    depths = sorted(reactions_by_depth.keys())

    # In retrosynthesis, we need to check if amine functionalization (later stage, lower depth)
    # is followed by nitro reduction (earlier stage, higher depth)
    for i in range(len(depths) - 1):
        later_stage_depth = depths[i]  # Lower depth number = later stage in synthesis
        earlier_stage_depth = depths[i + 1]  # Higher depth number = earlier stage in synthesis

        later_stage_reaction = reactions_by_depth[later_stage_depth]
        earlier_stage_reaction = reactions_by_depth[earlier_stage_depth]

        # Check if later stage reaction is amine functionalization
        # Common amine functionalization reactions include sulfonamide formation, amide formation, etc.
        amine_functionalization = False
        for rxn_type in AMINE_FUNCTIONALIZATION_REACTIONS:
            if checker.check_reaction(rxn_type, later_stage_reaction["rsmi"]):
                print(
                    f"Found amine functionalization reaction: {rxn_type} at depth {later_stage_depth}"
                )
                amine_functionalization = True
                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                break

        # Check if earlier stage reaction is nitro reduction
        nitro_reduction = False
        nitro_reduction_reaction_name = "Reduction of nitro groups to amines"
        if checker.check_reaction(
            nitro_reduction_reaction_name, earlier_stage_reaction["rsmi"]
        ):
            print(f"Found nitro reduction reaction at depth {earlier_stage_depth}")
            nitro_reduction = True
            if nitro_reduction_reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                findings_json["atomic_checks"]["named_reactions"].append(nitro_reduction_reaction_name)

        # If both conditions are met, we found the sequence
        if amine_functionalization and nitro_reduction:
            print(
                f"Found nitro reduction at depth {earlier_stage_depth} followed by amine functionalization at depth {later_stage_depth}"
            )
            result = True
            # Add the structural constraint to findings_json
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "Reduction of nitro groups to amines",
                    "after": {
                        "any_of": [
                            "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                            "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
                            "Acylation of primary amines",
                            "Acylation of secondary amines",
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N"
                        ]
                    },
                    "description": "Checks for a nitro reduction reaction that is followed in a subsequent step by an amine functionalization reaction."
                }
            })
            # Since we found the sequence, we can return early
            return result, findings_json

    return result, findings_json
