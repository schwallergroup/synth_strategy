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


# Refactored lists
REDUCTION_REACTIONS = [
    "Reduction of aldehydes and ketones to alcohols",
    "Reduction of carboxylic acid to primary alcohol",
    "Reduction of ester to primary alcohol",
    "Reduction of nitrile to amine",
    "Reduction of nitro groups to amines",
    "Reduction of primary amides to amines",
    "Reduction of secondary amides to amines",
    "Reduction of tertiary amides to amines",
    "Hydrogenation (double to single)",
    "Hydrogenation (triple to double)",
    "Arene hydrogenation",
    "Grignard from aldehyde to alcohol",
    "Grignard from ketone to alcohol",
]

OXIDATION_REACTIONS = [
    "Oxidation of aldehydes to carboxylic acids",
    "Oxidation of alkene to carboxylic acid",
    "Oxidation of alcohol to carboxylic acid",
    "Oxidation of ketone to carboxylic acid",
    "Oxidation of nitrile to carboxylic acid",
    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
    "Oxidation of alkene to aldehyde",
    "Oxidative esterification of primary alcohols",
    "Oxidation of alcohol and aldehyde to ester",
    "Quinone formation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Identifies synthesis routes containing at least one reduction and one oxidation step
    by matching each step against predefined lists of named redox reactions.
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

    reductions = 0
    oxidations = 0

    def dfs_traverse(node, depth=0):
        nonlocal reductions, oxidations, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for reduction reactions
            for reaction_type in REDUCTION_REACTIONS:
                if checker.check_reaction(reaction_type, rsmi):
                    reductions += 1
                    if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    break

            # Check for oxidation reactions
            for reaction_type in OXIDATION_REACTIONS:
                if checker.check_reaction(reaction_type, rsmi):
                    oxidations += 1
                    if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    break

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a chemical node
            next_depth = depth + 1

        # Traverse children with the determined next_depth
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least one reduction and one oxidation were found
    result = reductions >= 1 and oxidations >= 1

    if result:
        # Add the structural constraint if both reduction and oxidation are found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "any_reduction_reaction",
                    "any_oxidation_reaction"
                ],
                "min_counts": [
                    1,
                    1
                ]
            }
        })

    return result, findings_json
