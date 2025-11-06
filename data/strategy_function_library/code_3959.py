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
    This function detects if the synthetic route employs a late-stage nitro reduction strategy,
    where a nitro group is maintained through multiple steps and reduced to an amine in the final step.
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

    # Track if we've found a nitro reduction
    nitro_reduction_found = False
    # Track the depth at which nitro reduction occurs
    nitro_reduction_depth = None
    # Track if nitro group exists in earlier steps
    nitro_in_earlier_steps = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_found, nitro_reduction_depth, nitro_in_earlier_steps, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro group in reactants
            reactants_have_nitro = False
            for reactant in reactants_smiles:
                if reactant and checker.check_fg("Nitro group", reactant):
                    reactants_have_nitro = True
                    if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
                    break

            # Check for nitro group in product
            product_has_nitro = checker.check_fg("Nitro group", product_smiles)
            if product_has_nitro:
                if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

            # Check if this is a nitro reduction reaction
            is_nitro_reduction = checker.check_reaction("Reduction of nitro groups to amines", rsmi)

            if is_nitro_reduction:
                nitro_reduction_found = True
                nitro_reduction_depth = depth
                if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

            # Track if nitro exists in steps before the reduction
            if (reactants_have_nitro or product_has_nitro) and (
                nitro_reduction_depth is None or depth > nitro_reduction_depth
            ):
                nitro_in_earlier_steps = True

        elif node["type"] == "mol":
            # Check for nitro group in molecule nodes too
            if node["smiles"] and checker.check_fg("Nitro group", node["smiles"]):
                if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
                if nitro_reduction_depth is None or depth > nitro_reduction_depth:
                    nitro_in_earlier_steps = True

        # Recursively process children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # The strategy is present if we found a nitro reduction at a late stage (depth 0 or 1)
    # and the nitro group was present in earlier steps
    strategy_present = (
        nitro_reduction_found
        and nitro_in_earlier_steps
        and nitro_reduction_depth is not None
        and nitro_reduction_depth <= 1
    )

    if strategy_present:
        # Add structural constraints if the strategy is found
        # This corresponds to the 'sequence' constraint
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "before": {
                    "type": "functional_group",
                    "name": "Nitro group"
                },
                "after": {
                    "type": "named_reaction",
                    "name": "Reduction of nitro groups to amines"
                }
            }
        })
        # This corresponds to the 'positional' constraint
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Reduction of nitro groups to amines",
                "max_depth": 1
            }
        })

    return strategy_present, findings_json
