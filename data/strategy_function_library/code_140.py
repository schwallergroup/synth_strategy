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


OXIDATION_REACTIONS = [
    "Oxidation of aldehydes to carboxylic acids",
    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
    "Oxidation of alcohol to carboxylic acid",
    "Oxidation of ketone to carboxylic acid",
    "Oxidation of alkene to carboxylic acid",
    "Oxidation of nitrile to carboxylic acid",
    "Oxidation of amide to carboxylic acid",
    "Oxidation of alkene to aldehyde",
    "Oxidative esterification of primary alcohols",
    "Oxidation of alcohol and aldehyde to ester",
    "Quinone formation",
    "Aromatic hydroxylation",
    "Sulfanyl to sulfinyl_peroxide",
    "Sulfanyl to sulfinyl_H2O2",
    "Sulfanyl to sulfinyl",
    "Aerobic oxidation of Grignard reagents",
]

REDUCTION_REACTIONS = [
    "Azide to amine reduction (Staudinger)",
    "Reduction of aldehydes and ketones to alcohols",
    "Reduction of ester to primary alcohol",
    "Reduction of ketone to secondary alcohol",
    "Reduction of carboxylic acid to primary alcohol",
    "Reduction of nitro groups to amines",
    "Reduction of primary amides to amines",
    "Reduction of secondary amides to amines",
    "Reduction of tertiary amides to amines",
    "Reduction of nitrile to amine",
    "Hydrogenation (double to single)",
    "Hydrogenation (triple to double)",
    "Arene hydrogenation",
    "Reductive amination with aldehyde",
    "Reductive amination with ketone",
    "Reductive amination with alcohol",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if a synthesis route contains at least one reaction from a defined set of oxidation transformations AND at least one reaction from a defined set of reduction transformations.
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

    # Track oxidation and reduction steps
    found_oxidation = False
    found_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_oxidation, found_reduction, findings_json

        # Terminate early if both are already found
        if found_oxidation and found_reduction:
            return

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for oxidation reactions if not already found
            if not found_oxidation:
                for r in OXIDATION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        found_oxidation = True
                        findings_json["atomic_checks"]["named_reactions"].append(r)
                        break

            # Check for reduction reactions if not already found
            if not found_reduction:
                for r in REDUCTION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        found_reduction = True
                        findings_json["atomic_checks"]["named_reactions"].append(r)
                        break

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # This means it's a 'chemical' node
            next_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            # Stop traversing this branch if both have been found
            if found_oxidation and found_reduction:
                break
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Return True if both oxidation and reduction steps are found
    result = found_oxidation and found_reduction

    if result:
        # Add the structural constraint if both oxidation and reduction were found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "any_oxidation_reaction",
                    "any_reduction_reaction"
                ]
            }
        })

    return result, findings_json
